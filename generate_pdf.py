#!/usr/bin/env python3
"""Generate PDF white paper with embedded figures from markdown."""

import os
import re
import markdown
from weasyprint import HTML, CSS

os.chdir('/home/user/abc_atlas_access')

# Read the markdown
with open('white_paper.md', 'r') as f:
    md_content = f.read()

# Convert markdown to HTML
html_body = markdown.markdown(md_content, extensions=['tables', 'toc'])

# Fix image paths to be absolute for weasyprint
html_body = html_body.replace('src="outputs/', f'src="file:///home/user/abc_atlas_access/outputs/')

# Wrap in full HTML with styling
html_doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
@page {{
    size: letter;
    margin: 1in 0.75in;
    @bottom-center {{
        content: counter(page);
        font-size: 10pt;
        color: #666;
    }}
}}
body {{
    font-family: "Times New Roman", Times, Georgia, serif;
    font-size: 11pt;
    line-height: 1.5;
    color: #000;
    max-width: 100%;
}}
h1 {{
    font-size: 16pt;
    text-align: center;
    margin-bottom: 0.5em;
    line-height: 1.3;
    page-break-after: avoid;
}}
h2 {{
    font-size: 14pt;
    margin-top: 1.5em;
    margin-bottom: 0.5em;
    page-break-after: avoid;
}}
h3 {{
    font-size: 12pt;
    margin-top: 1em;
    margin-bottom: 0.5em;
    page-break-after: avoid;
}}
p {{
    text-align: justify;
    margin-bottom: 0.8em;
}}
img {{
    max-width: 100%;
    height: auto;
    display: block;
    margin: 1em auto;
    page-break-inside: avoid;
}}
hr {{
    border: none;
    border-top: 1px solid #ccc;
    margin: 1.5em 0;
}}
strong {{
    font-weight: bold;
}}
em {{
    font-style: italic;
}}
</style>
</head>
<body>
{html_body}
</body>
</html>
"""

# Write HTML for debugging
with open('white_paper.html', 'w') as f:
    f.write(html_doc)

print("Converting to PDF...")
HTML(string=html_doc, base_url='/home/user/abc_atlas_access/').write_pdf(
    'white_paper.pdf',
)
print("PDF generated: white_paper.pdf")
print(f"Size: {os.path.getsize('white_paper.pdf') / 1024 / 1024:.1f} MB")
