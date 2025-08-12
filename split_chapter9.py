#!/usr/bin/env python3
"""
Script to split Chapter 9 (Stream Restoration) into manageable chunks for processing.
Based on the successful chunking approach used for previous chapters.
"""

import fitz  # PyMuPDF
import os
from pathlib import Path

def split_chapter9(input_pdf_path, output_dir, pages_per_chunk=5):
    """Split Chapter 9 PDF into smaller chunks of specified page count."""
    
    if not os.path.exists(input_pdf_path):
        print(f"PDF file not found: {input_pdf_path}")
        return
    
    # Open the PDF
    doc = fitz.open(input_pdf_path)
    total_pages = len(doc)
    
    print(f"Chapter 9 has {total_pages} pages")
    print(f"Splitting into chunks of {pages_per_chunk} pages each")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    chunk_num = 1
    start_page = 0
    
    while start_page < total_pages:
        # Calculate end page for this chunk
        end_page = min(start_page + pages_per_chunk, total_pages)
        
        # Create new PDF for this chunk
        chunk_doc = fitz.open()
        
        # Add pages to chunk
        for page_num in range(start_page, end_page):
            chunk_doc.insert_pdf(doc, from_page=page_num, to_page=page_num)
        
        # Save chunk
        chunk_filename = f"chapter09_chunk_{chunk_num:02d}_pages_{start_page+1}-{end_page}.pdf"
        chunk_path = output_path / chunk_filename
        chunk_doc.save(str(chunk_path))
        chunk_doc.close()
        
        print(f"  Chunk {chunk_num}: Pages {start_page+1}-{end_page} -> {chunk_filename}")
        
        # Move to next chunk
        start_page = end_page
        chunk_num += 1
    
    doc.close()
    
    print(f"\nSuccessfully split Chapter 9 into {chunk_num-1} chunks")
    return chunk_num - 1

def main():
    chapter9_pdf = "chapters/Chapter_12_CHAPTER_9_STREAM_RESTORATION.pdf"
    output_directory = "chapters/chapter09_chunks"
    
    print("Splitting Chapter 9 (Stream Restoration)...")
    
    # Split into 5-page chunks
    num_chunks = split_chapter9(chapter9_pdf, output_directory, pages_per_chunk=5)
    
    if num_chunks:
        print(f"\nChapter 9 successfully split into {num_chunks} chunks in {output_directory}/")
        print("Ready for sequential analysis and implementation planning.")

if __name__ == "__main__":
    main()