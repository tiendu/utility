import requests
import re
import os
import time

# Base URL for vector search
base_url = 'https://www.addgene.org/search/catalog/plasmids/?q=&page_size=50&vector_types=AAV&page_number='
page_number = 1

# Ensure the output directory exists
output_dir = 'plasmid_sequences'
os.makedirs(output_dir, exist_ok=True)

# Function to get total number of results
def get_total_results():
    initial_response = requests.get(base_url + str(1))
    total_results_match = re.search(r'<p class="results-text">([0-9,]+) results</p>', initial_response.text)
    if total_results_match:
        return int(total_results_match.group(1).replace(',', ''))
    return 0

total_results = get_total_results()
collected_plasmid_ids = set()

while len(collected_plasmid_ids) < total_results:
    search_url = base_url + str(page_number)
    search_response = requests.get(search_url)

    if search_response.status_code != 200:
        print("Failed to fetch:", search_url)
        break

    # Extract plasmid IDs from the search results page
    plasmid_ids = re.findall(r'<div id="Plasmids-([0-9]+)" class="search-result-item">', search_response.text)
    for plasmid_id in plasmid_ids:
        if plasmid_id in collected_plasmid_ids:
            continue
        collected_plasmid_ids.add(plasmid_id)
        sequence_url = f'https://www.addgene.org/{plasmid_id}/sequences/'
        sequence_response = requests.get(sequence_url)

        if sequence_response.status_code != 200:
            print(f"Failed to fetch sequence page for plasmid ID: {plasmid_id}")
            continue

        # Find the .gbk download link
        gbk_link = re.search(r'<a href="(https://media\.addgene\.org/[^"]*\.gbk)"', sequence_response.text)
        if gbk_link:
            gbk_url = gbk_link.group(1)
            gbk_filename = os.path.join(output_dir, f"plasmid_{plasmid_id}.gbk")
            gbk_content = requests.get(gbk_url).content

            with open(gbk_filename, 'wb') as f:
                f.write(gbk_content)
            print(f"Downloaded GBK file: {gbk_filename}")
        else:
            print(f"No GBK file found for plasmid ID: {plasmid_id}")

    # Check if there is a next page, and reset if needed
    if '<li id="next-btn" class="page-item disabled">' in search_response.text:
        if page_number == 1:
            break
        else:
            page_number = 1  # Restart from the first page to catch any shuffled sequences
            time.sleep(2)  # Delay to avoid hammering the server
    else:
        page_number += 1

print(f"Scraping and downloading complete. Total unique plasmid IDs collected: {len(collected_plasmid_ids)}")