import requests
from bs4 import BeautifulSoup
import html2text
import os
import time

BASE_URL = "https://bookdown.org/sjcockell/ismb-tutorial-2023/"
PAGES = [
    "index.html",
    "practical-session-1.html",
    "practical-session-2.html",
    "practical-session-3.html",
    "practical-session-4.html",
    "references.html"
]

def scrape_page(page_name):
    url = f"{BASE_URL}{page_name}"
    print(f"Scraping {url}...")
    try:
        response = requests.get(url)
        response.raise_for_status()
    except Exception as e:
        print(f"Error fetching {url}: {e}")
        return ""

    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Bookdown main content is typically in a <section class="normal"> or .page-inner
    content = soup.find('section', class_='normal')
    if not content:
        content = soup.find('div', class_='page-inner')
    
    if not content:
        print(f"Could not find main content for {page_name}")
        return ""

    # Remove navigation arrows if present
    for nav in content.find_all('a', class_='navigation'):
        nav.decompose()

    h = html2text.HTML2Text()
    h.ignore_links = False
    h.body_width = 0  # No wrapping
    
    markdown = h.handle(str(content))
    return markdown

def main():
    full_markdown = ""
    for page in PAGES:
        markdown = scrape_page(page)
        if markdown:
            full_markdown += f"\n\n<!-- PAGE: {page} -->\n\n"
            full_markdown += markdown
        time.sleep(1)  # Respectful delay

    with open("scraped_course.md", "w") as f:
        f.write(full_markdown)
    print("Scraping complete. Saved to scraped_course.md")

if __name__ == "__main__":
    main()
