import requests
import json
import re
import sys
import csv
from bs4 import BeautifulSoup

def get_orpha_data(orpha_id):
    url = f"https://www.orpha.net/en/disease/detail/{orpha_id}"
    response = requests.get(url)
    if response.status_code == 200:
        html_content = response.text
        soup = BeautifulSoup(html_content, 'html.parser')

        # Extract the disease name from the title
        disease_name = soup.title.string.replace("Orphanet: ", "").strip().replace(" ", "%20")

        # Use regex to find the JSON-LD script
        json_ld_script = re.search(r'<script type="application/ld\+json">\s*(\{.*?\})\s*</script>', html_content, re.DOTALL)
        if json_ld_script:
            json_data = json_ld_script.group(1)
            data = json.loads(json_data)
            # Extract orphan drugs information using regex
            drugs_info = re.search(
                r'<a href="/en/drug\?orphaCode=\d+&amp;diseaseName=.*?&amp;mode=all">\s*Orphan designation\(s\) and orphan drug\(s\) \((\d+)\)\s*</a>',
                html_content
            )
            orphan_drugs_count = drugs_info.group(1) if drugs_info else "NA"

            orphan_drug_links = []
            if drugs_info:
                drugs_url = f"https://www.orpha.net/en/drug?orphaCode={orpha_id}&diseaseName={disease_name}&mode=all"
                orphan_drug_links = get_orphan_drug_links(drugs_url)

            return data, disease_name.replace("%20", " "), orphan_drugs_count, orphan_drug_links
        else:
            print(f"No JSON-LD script found for ID {orpha_id}")
            return None, disease_name.replace("%20", " "), "NA", []
    else:
        print(f"Failed to retrieve page for ID {orpha_id}")
        return None, "NA", "NA", []

def get_orphan_drug_links(drugs_url):
    response = requests.get(drugs_url)
    orphan_drug_links = []
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        drug_cards = soup.find_all("div", class_="drug-card")
        for card in drug_cards:
            drug_link = card.find("a", href=re.compile(r'/en/drug/substance/'))
            if drug_link:
                link = f"https://www.orpha.net{drug_link['href']}"
                orphan_drug_links.append(link)
    return orphan_drug_links

def get_drug_details(drug_url):
    response = requests.get(drug_url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Finding the drug name from h2 tag in div with class "result-detail"
        result_detail = soup.find("div", class_="result-detail")
        drug_name = result_detail.find("h2").text.strip() if result_detail else "NA"
        
        # Finding code/synonyms and chemical description using more precise searching
        code_synonyms = "NA"
        chemical_description = "NA"
        
        col6_tags = soup.find_all("div", class_="col-6")
        for tag in col6_tags:
            p_tags = tag.find_all("p")
            for p_tag in p_tags:
                if "Code/Synonyms:" in p_tag.text:
                    code_synonyms = p_tag.find("span", class_="fw-bold").text.strip()
                if "Chemical name or description:" in p_tag.text:
                    chemical_description = p_tag.find("span", class_="fw-bold").text.strip()

        return drug_name, code_synonyms, chemical_description
    return None, None, None

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

def main(file_path, output_csv):
    # Read Orpha IDs from a text file
    orpha_ids = read_orpha_ids(file_path)

    # Open CSV file for writing
    with open(output_csv, mode='w', newline='') as csv_file:
        fieldnames = ["Orpha ID", "GARD ID", "Disease Name", "USA Estimate", "World Estimate", "Orphan Drugs", "Drug Name", "Code/Synonyms", "Chemical Description"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        # Iterate over each Orpha ID and scrape the data
        for orpha_id in orpha_ids:
            orpha_data, disease_name, orphan_drugs_count, orphan_drug_links = get_orpha_data(orpha_id)
            if orpha_data:
                # Extract the GARD codeValue
                gard_code_value = "NA"
                for code in orpha_data.get("code", []):
                    if code.get("codingSystem") == "GARD":
                        gard_code_value = code.get("codeValue")
                        break
                
                usa_estimate = "NA"
                world_estimate = "NA"

                if gard_code_value != "NA":
                    gard_data = get_gard_data(gard_code_value)
                    if gard_data:
                        usa_estimate = gard_data.get("USA_Estimate__c", "NA")
                        world_estimate = gard_data.get("World_Estimate__c", "NA")
                
                # Truncate the disease name if it's too long
                disease_name_truncated = truncate_string(disease_name, 50)
                
                # Get details for each orphan drug
                for drug_url in orphan_drug_links:
                    drug_name, code_synonyms, chemical_description = get_drug_details(drug_url)
                    if drug_name:
                        writer.writerow({
                            "Orpha ID": orpha_id,
                            "GARD ID": gard_code_value,
                            "Disease Name": disease_name_truncated,
                            "USA Estimate": usa_estimate,
                            "World Estimate": world_estimate,
                            "Orphan Drugs": orphan_drugs_count,
                            "Drug Name": drug_name,
                            "Code/Synonyms": code_synonyms,
                            "Chemical Description": chemical_description
                        })
            else:
                # Append default NA values for missing orpha data
                writer.writerow({
                    "Orpha ID": orpha_id,
                    "GARD ID": "NA",
                    "Disease Name": disease_name,
                    "USA Estimate": "NA",
                    "World Estimate": "NA",
                    "Orphan Drugs": "NA",
                    "Drug Name": "NA",
                    "Code/Synonyms": "NA",
                    "Chemical Description": "NA"
                })

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <file_path> <output_csv>")
        sys.exit(1)
    file_path = sys.argv[1]
    output_csv = sys.argv[2]
    main(file_path, output_csv)
