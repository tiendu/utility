import requests
import json
import re
import sys
from tabulate import tabulate

def get_orpha_data(orpha_id):
    url = f"https://www.orpha.net/en/disease/detail/{orpha_id}"
    response = requests.get(url)
    if response.status_code == 200:
        html_content = response.text
        # Use regex to find the JSON-LD script
        json_ld_script = re.search(r'<script type="application/ld\+json">\s*(\{.*?\})\s*</script>', html_content, re.DOTALL)
        if json_ld_script:
            json_data = json_ld_script.group(1)
            data = json.loads(json_data)
            return data
        else:
            print(f"No JSON-LD script found for ID {orpha_id}")
            return None
    else:
        print(f"Failed to retrieve page for ID {orpha_id}")
        return None

def get_gard_data(code_value):
    url = f"https://rarediseases.info.nih.gov/assets/singles/{code_value}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Failed to retrieve GARD data for code {code_value}")
        return None

def read_orpha_ids(file_path):
    with open(file_path, 'r') as file:
        orpha_ids = file.read().splitlines()
    return orpha_ids

def truncate_string(string: str, max_length: int) -> str:
    if len(string) > max_length:
        return string[:max_length - 3] + '...'
    return string

def main(file_path):
    # Read Orpha IDs from a text file
    orpha_ids = read_orpha_ids(file_path)

    # Initialize an empty list to store the results
    results = []

    # Iterate over each Orpha ID and scrape the data
    for orpha_id in orpha_ids:
        orpha_data = get_orpha_data(orpha_id)
        if orpha_data:
            # Extract the GARD codeValue
            gard_code_value = "NA"
            for code in orpha_data.get("code", []):
                if code.get("codingSystem") == "GARD":
                    gard_code_value = code.get("codeValue")
                    break
            
            usa_estimate = "NA"
            world_estimate = "NA"
            disease_name = "NA"

            if gard_code_value != "NA":
                gard_data = get_gard_data(gard_code_value)
                if gard_data:
                    usa_estimate = gard_data.get("USA_Estimate__c", "NA")
                    world_estimate = gard_data.get("World_Estimate__c", "NA")
                    disease_name = gard_data.get("Disease_Name_Full__c", "NA")
            
            # Truncate the disease name if it's too long
            disease_name_truncated = truncate_string(disease_name, 50)
            
            # Append the result to the list
            results.append({
                "Orpha ID": orpha_id,
                "GARD ID": gard_code_value,
                "Disease Name": disease_name_truncated,
                "USA Estimate": usa_estimate,
                "World Estimate": world_estimate
            })
        else:
            # Append default NA values for missing orpha data
            results.append({
                "Orpha ID": orpha_id,
                "GARD ID": "NA",
                "Disease Name": "NA",
                "USA Estimate": "NA",
                "World Estimate": "NA"
            })

    # Create a table from the results
    table = [
        ["Orpha ID", "GARD ID", "Disease Name", "USA Estimate", "World Estimate"]
    ]

    for result in results:
        table.append([result["Orpha ID"], result["GARD ID"], result["Disease Name"], result["USA Estimate"], result["World Estimate"]])

    # Display the table using tabulate
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)
    file_path = sys.argv[1]
    main(file_path)
