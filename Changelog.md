# Changelog

All notable changes to the **Lemos Lambda** project will be documented in this file.

## [1.6.1] - 2025-12-19
### Added (Architecture of Silence)
- **Functional Context Filter:** Implemented a pre-screening layer that rejects articles lacking mechanistic vocabulary (e.g., "signaling", "pathway", "relaxation"), effectively filtering out purely clinical, epidemiological, or administrative records.
- **Deep Mining Engine:** Replaced the lightweight scraper with a robust Biopython-based engine capable of analyzing batch abstracts (n=60) with relevance sorting.
- **Semantic Wall:** Introduced a multi-layer filtration system:
    - **Layer 1 (Query):** Restricted to `[Title/Abstract]` fields to avoid affiliation noise.
    - **Layer 2 (Exclusion):** Explicit boolean NOT for Reviews and Editorials.
    - **Layer 3 (Regex):** Strict alphanumeric enforcement (requires digits or hyphens for validity).
    - **Layer 4 (Blacklist):** Added aggressive exclusion for radiotherapy acronyms (PTV, GTV), grant codes (DK...), and administrative metadata (EDAT, MHDA).
- **Inverted Lambda Score:** Mathematics flipped to highlight "Blue Ocean" opportunities. The score now represents (Global Popularity / Local Saturation), prioritizing high-impact/low-competition targets.

### Changed
- **Blue Ocean Button:** Now dynamic. Instead of loading a static list, it triggers a real-time mining session on PubMed based on the user's input organ.
- **Auto-Detect Logic:** Refined to use the "Semantic Wall" logic, reducing false positives (e.g., stopping "PUMCH" or "OABSS" from being flagged as genes).
- **UI Adjustments:** Removed donation links to maintain academic neutrality.

### Fixed
- **Grant Code Noise:** Fixed regex leakage where grant numbers (e.g., KMUH-110) were being detected as protein isoforms.
- **Metadata Leak:** Fixed issue where PubMed XML tags (e.g., PHST, AID) appeared in results.