# Non-Rigid-Shape-Correspondence

**Unsupervised Landmark Discovery via Karcher Means for Non-Rigid Shape Correspondence**
**Sample**: Ordering Landmarks on Non-Uniform Mesh Surfaces
![demo image of computing ordered landmarks on non-uniform mesh surfaces ](https://github.com/ZHLRJ/Non-Rigid-Shape-Correspondence/blob/main/teaser_cloud.png?raw=true)


This repository contains the official implementation of our method for discovering intrinsic landmarks in 3D shapes for dense correspondence tasks. Our approach computes unsupervised, repeatable, and geometrically meaningful landmarks based on the intrinsic geometry of the surface, enabling robust shape matching without the need for manually labeled annotations.

## 🔍 Overview

Establishing dense shape correspondence across non-rigid 3D surfaces is a long-standing challenge. This project introduces an unsupervised method that:
- Computes intrinsic landmarks using **geometric medians** (Karcher means),
- Selects salient points via **farthest point sampling** over the intrinsic domain,
- Supports matching under **near-isometric deformations**,
- Can be integrated into existing **functional map pipelines**.


## 📁 Project Structure
```text
Non-Rigid-Shape-Correspondence/
├── geometricMedian/     # Karcher mean computation for intrinsic landmark localization
├── logarithmicMap/      # Intrinsic logarithmic map computation on triangle meshes
├── pyFM/                # Functional map pipeline and shape descriptor tools
├── tracinggeodesic/     # Geodesic tracing and shortest path tools
└── README.md            # Project documentation and usage guide

