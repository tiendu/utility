import requests
import re

# Base URL for vector search
base_url = 'https://www.addgene.org/search/catalog/plasmids/?q=&page_size=50&vector_types=AAV&page_number='
page_number = 1

while True:
    search_url = base_url + str(page_number)
    search_response = requests.get(search_url)

    plasmid_ids = re.findall(r'<div id="Plasmids-([0-9]+)" class="search-result-item">', search_response.text)
    for plasmid_id in plasmid_ids:
        sequence_url = 'https://www.addgene.org/' + str(plasmid_id) + '/sequences/'
        sequence_response = requests.get(sequence_url)
        gbk_link = re.search(r'<a href="(https://media\.addgene\.org/[^"]*\.gbk)"', sequence_response.text)
        if gbk_link:
            gbk_url = gbk_link.group(1)
            gbk_filename = f"plasmid_{plasmid_id}.gbk"
            gbk_content = requests.get(gbk_url).content
            with open(gbk_filename, 'wb') as f:
                f.write(gbk_content)
            print(f"Downloaded GBK file: {gbk_filename}")

    if search_response.status_code == 200:
        if '<li id="next-btn" class="page-item disabled">' in search_response.text:
            break
        else:
            page_number += 1
    else:
        print("Failed to fetch:", search_url)
        break
