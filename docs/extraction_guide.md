# PDF Processing and Verification Workflow

This document outlines the standard operating procedure for generating or verifying course materials from PDF lectures using image extraction.

## Goal
To synthesize accurate Markdown study guides or verify notes by processing PDF slides as high-fidelity images.

## Procedure

### 1. Preparation
1.  **Identify Target Files**: Locate source PDF(s) in `_pdf/`.
2.  **Create Image Directory**: Create a temporary directory (e.g., `_temp_images/`) for processing.

### 2. Image Extraction
Convert PDF pages to JPEG images for direct visual analysis by the agent. This bypasses text extraction issues and layout distortions.

**Execution:**
```bash
pdftoppm -jpeg -r 150 source_lecture.pdf _temp_images/page
```
*Note: `-r 150` provides sufficient resolution while keeping file sizes manageable.*

### 3. Iterative Batch Analysis & Implementation
**CRITICAL**: Process images in batches of 5. Update target file immediately after each batch to prevent context loss.

**Step-by-Step Loop:**
1.  **Batch Load**: Read **5** consecutive page images (e.g., `page-01.jpg` to `page-05.jpg`) using `read_file`.
2.  **Analyze & Commit**:
    -   **Extract Content**: Identify headers, key concepts, formulas, and code snippets.
    -   **Write/Update**: Immediately append synthesized content to the target Markdown file or update existing sections. Do not wait for all batches to finish.
    -   **Cross-Reference (if verifying)**: Compare image content against existing Markdown and note discrepancies.
3.  **Repeat**: Move to the next batch of 5 images until the lecture is complete.

### 4. Final Review
-   Review the file for cohesive flow and check for missing transitions between batches.
-   Ensure consistent formatting throughout.
-   Add WikiLinks to related materials.

### 5. Cleanup
-   Remove the `_temp_images/` directory when finished.

---
**Instruction to Agent:**
Follow this procedure for study guide generation or verification. Use `pdftoppm` to extract images and process them in batches of 5. Do not use split PDF scripts or Ghostscript.