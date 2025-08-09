#!/usr/bin/env python3
"""Split Chapter 4 PDF into 5-page chunks."""

import fitz
import math
import os

def split_pdf_into_chunks(input_path, output_dir, pages_per_chunk=5):
    """Split a PDF into smaller chunks of specified page count."""
    
    # Open the PDF
    doc = fitz.open(input_path)
    total_pages = len(doc)
    
    print(f"Total pages in Chapter 4: {total_pages}")
    
    # Calculate number of chunks needed
    num_chunks = math.ceil(total_pages / pages_per_chunk)
    print(f"Creating {num_chunks} chunks of {pages_per_chunk} pages each")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Split into chunks
    for chunk_num in range(num_chunks):
        start_page = chunk_num * pages_per_chunk
        end_page = min(start_page + pages_per_chunk - 1, total_pages - 1)
        
        # Create a new document for this chunk
        chunk_doc = fitz.open()
        
        # Insert pages from original document
        chunk_doc.insert_pdf(doc, from_page=start_page, to_page=end_page)
        
        # Save the chunk
        chunk_filename = f"chapter04_chunk_{chunk_num + 1:02d}_pages_{start_page + 1}-{end_page + 1}.pdf"
        chunk_path = os.path.join(output_dir, chunk_filename)
        chunk_doc.save(chunk_path)
        chunk_doc.close()
        
        print(f"Created: {chunk_filename} (pages {start_page + 1}-{end_page + 1})")
    
    doc.close()
    print(f"Successfully split Chapter 4 into {num_chunks} chunks")

if __name__ == "__main__":
    input_pdf = "chapters/Chapter_07_CHAPTER_4_FINE_GRAINED_SEDIMENT_TRANSPORT.pdf"
    output_directory = "chapters/chapter04_chunks"
    
    split_pdf_into_chunks(input_pdf, output_directory, pages_per_chunk=5)