import os
import zipfile
import argparse
import getpass
import subprocess

def install_required_packages():
    required_packages = ["PyGithub"]
    installed = subprocess.check_output(['pip', 'freeze']).decode().splitlines()
    installed_packages = [pkg.split('==')[0] for pkg in installed]
    for package in required_packages:
        if package not in installed_packages:
            subprocess.check_call(["pip", "install", package])

install_required_packages()

from github import Github

def backup_github_to_single_zip(token, output_zip='all_repos.zip'):
    # GitHub setup
    g = Github(token)
    user = g.get_user()

    # List all your repos
    repos = user.get_repos()

    for repo in repos:
        repo_name = repo.name
        print(f"Cloning {repo_name}...")

        # Clone the repo using the token for authentication
        modified_clone_url = repo.clone_url.replace('https://', f'https://{token}@')
        os.system(f"git clone {modified_clone_url}")

    # Zip all the repo folders into a single zip file
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for repo in repos:
            for foldername, subfolders, filenames in os.walk(repo.name):
                for filename in filenames:
                    zipf.write(os.path.join(foldername, filename))

    print(f"Backup completed! All repos zipped into: {output_zip}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Backup GitHub repositories into a single ZIP file.")
    parser.add_argument("--token", type=str, help="Your GitHub Personal Access Token.")
    args = parser.parse_args()

    if not args.token:
        args.token = getpass.getpass(prompt='Enter your GitHub token (input will be hidden): ')

    backup_github_to_single_zip(args.token)
