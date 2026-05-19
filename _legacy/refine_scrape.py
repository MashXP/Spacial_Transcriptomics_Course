import requests
from bs4 import BeautifulSoup
import html2text
import os
import time
import re

BASE_URL = "https://bookdown.org/sjcockell/ismb-tutorial-2023/"
PAGES = [
    "index.html",
    "practical-session-1.html",
    "practical-session-2.html",
    "practical-session-3.html",
    "practical-session-4.html",
    "references.html"
]

def clean_markdown(text):
    # Remove empty brackets like [] used for anchor links
    text = re.sub(r'\[\]\([^)]*\)', '', text)
    # Remove solitary []
    text = re.sub(r'\[\]', '', text)
    # Fix weird character issues (e.g., â)
    text = text.replace('â', '"').replace('âs', "'s").replace('â', "'")
    # Convert [code] ... [/code] to ```R ... ```
    # Using non-greedy match to handle multiple blocks on one page
    text = re.sub(r'\[code\]\s*', '```R\n', text)
    text = re.sub(r'\s*\[/code\]', '\n```', text)
    # Ensure code blocks have proper spacing
    text = re.sub(r'```\n\n', '```\n', text)
    text = re.sub(r'\n\n```', '\n```', text)
    # Remove leading spaces from code block lines (often 4 spaces in bookdown output)
    def fix_code_indent(match):
        block = match.group(1)
        lines = block.split('\n')
        fixed_lines = [line.lstrip() for line in lines]
        return '```R\n' + '\n'.join(fixed_lines) + '\n```'
    
    text = re.compile(r'```R\n(.*?)\n```', re.DOTALL).sub(fix_code_indent, text)

    return text

def scrape_and_process():
    if not os.path.exists("course_4_split"):
        os.makedirs("course_4_split")

    for page in PAGES:
        url = f"{BASE_URL}{page}"
        print(f"Processing {url}...")
        try:
            response = requests.get(url)
            response.raise_for_status()
            # Fix encoding
            response.encoding = response.apparent_encoding
        except Exception as e:
            print(f"Error fetching {url}: {e}")
            continue

        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Bookdown content
        content = soup.find('section', class_='normal')
        if not content:
            content = soup.find('div', class_='page-inner')
        
        if not content:
            print(f"Could not find content for {page}")
            continue

        # Remove nav
        for nav in content.find_all('a', class_='navigation'):
            nav.decompose()
        
        # Remove the 'header-section-number' spans if they clutter
        for span in content.find_all('span', class_='header-section-number'):
            # Keep the text but maybe wrap it? Or just leave it. 
            # Often they look like "1.1"
            pass

        h = html2text.HTML2Text()
        h.ignore_links = False
        h.body_width = 0
        h.mark_code = True # Try to improve code block detection
        
        raw_md = h.handle(str(content))
        cleaned_md = clean_markdown(raw_md)
        
        file_name = page.replace(".html", ".md")
        output_path = os.path.join("course_4_split", file_name)
        
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(cleaned_md)
        
        print(f"Saved {output_path}")
        time.sleep(0.5)

if __name__ == "__main__":
    scrape_and_process()
