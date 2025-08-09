#!/usr/bin/env python3
"""
Script to analyze and split the large sedimentation engineering PDF by chapters.
"""

import fitz  # PyMuPDF
import os
import re
from pathlib import Path

def analyze_pdf_structure(pdf_path):
    """Analyze PDF to find chapter boundaries."""
    doc = fitz.open(pdf_path)
    
    print(f"PDF has {len(doc)} pages")
    
    # Get table of contents
    toc = doc.get_toc()
    if toc:
        print("\nTable of Contents found:")
        for level, title, page in toc[:20]:  # Show first 20 entries
            indent = "  " * (level - 1)
            print(f"{indent}{page}: {title}")
    
    # Look for chapter patterns in first few pages
    chapter_pattern = re.compile(r'(?i)chapter\s+\d+|part\s+\d+|section\s+\d+', re.MULTILINE)
    
    chapters = []
    for page_num in range(min(50, len(doc))):  # Check first 50 pages for structure
        page = doc[page_num]
        text = page.get_text()
        
        # Look for chapter headings
        for match in chapter_pattern.finditer(text):
            context = text[max(0, match.start()-50):match.end()+50].strip()
            if len(context) < 200:  # Likely a heading, not body text
                chapters.append((page_num + 1, match.group(), context))
    
    print(f"\nFound {len(chapters)} potential chapter markers:")
    for page, title, context in chapters[:10]:  # Show first 10
        print(f"Page {page}: {title}")
        print(f"  Context: {context[:100]}...")
    
    doc.close()
    return toc, chapters

def split_pdf_by_toc(pdf_path, output_dir):
    """Split PDF using table of contents."""
    doc = fitz.open(pdf_path)
    toc = doc.get_toc()
    
    if not toc:
        print("No table of contents found")
        return
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Filter for main chapters (level 1 entries)
    main_chapters = [(title, page) for level, title, page in toc if level == 1]
    
    print(f"Splitting into {len(main_chapters)} main chapters:")
    
    for i, (title, start_page) in enumerate(main_chapters):
        # Clean title for filename
        clean_title = re.sub(r'[^\w\s-]', '', title).strip()
        clean_title = re.sub(r'[-\s]+', '_', clean_title)
        
        # Determine end page
        if i + 1 < len(main_chapters):
            end_page = main_chapters[i + 1][1] - 1
        else:
            end_page = len(doc)
        
        # Create new PDF for this chapter
        chapter_doc = fitz.open()
        for page_num in range(start_page - 1, min(end_page, len(doc))):
            chapter_doc.insert_pdf(doc, from_page=page_num, to_page=page_num)
        
        output_path = output_dir / f"Chapter_{i+1:02d}_{clean_title}.pdf"
        chapter_doc.save(str(output_path))
        chapter_doc.close()
        
        print(f"  Chapter {i+1}: {title} (Pages {start_page}-{end_page}) -> {output_path.name}")
    
    doc.close()

def main():
    pdf_path = "Environmental and Water Resources Institute (U.S.) and García - 2008 - Sedimentation engineering processes, measurements.pdf"
    
    if not os.path.exists(pdf_path):
        print(f"PDF file not found: {pdf_path}")
        return
    
    print("Analyzing PDF structure...")
    toc, chapters = analyze_pdf_structure(pdf_path)
    
    if toc:
        print("\nSplitting PDF by table of contents...")
        split_pdf_by_toc(pdf_path, "chapters")
    else:
        print("\nNo table of contents found. Manual chapter detection would be needed.")

if __name__ == "__main__":
    main()