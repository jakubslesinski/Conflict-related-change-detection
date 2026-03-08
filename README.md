# 📄 Integrating optical and radar satellite data for conflict-related change detection in Ukraine: a multi-temporal analysis of building destruction and agricultural disruption

---

<a href="./figs/Background_pwtt.png" target="_blank">
  <img width="100%" src="./figs/Background_pwtt.png" alt="Background image">
</a>

---

## 📝 Abstract
The ongoing war in Ukraine has caused extensive damage to infrastructure, agriculture and the environment, while ground-based evaluation remains severely limited due to safety concerns. This paper introduces a novel methodology for change detection based solely on open-access Sentinel-1 and Sentinel-2 data. The approach integrates multi-temporal optical image classification, enhanced by a conditional correction mechanism, with radar-based change detection methods (REACTIV, PWTT, Omnibus). In addition to annual assessments, the method incorporates seasonal analyses of agricultural activity, enabling a more detailed evaluation of war-related disruptions. This fusion allowed for a multi-faceted assessment of building destruction and agricultural decline between 2022 and 2025, with particular emphasis on the Bakhmut region. A comparison of the results of the proposed method with the UNOSAT database confirmed more than 80% of building damage. The proposed classification methodology achieved higher precision in built-up area detection compared to the AlphaEarth platform (0.98 vs. 0.87). The results demonstrate that combining optical and radar techniques provides a robust and scalable tool for monitoring the impacts of armed conflicts. The code are available at:
https://github.com/jakubslesinski/Conflict-related-change-detection

## 📖 Citation (to be completed) 

If you use this repository, please cite:

> Surname, N., & Coauthor, M. (2025). *Article Title*. Journal Name. https://doi.org/xxxx  

```bibtex
@article{surname2025article,
  title   = {Integrating optical and radar satellite data for conflict-related change detection in Ukraine},
  author  = {Kinga Karwowska, Jakub Slesinski, Aleksandra Sekrecka, Michal Smiarowski, Kart Metsoja},
  journal = {Scientific Reports},
  year    = {2026},
  doi     = {https://doi.org/10.1038/s41598-026-41424-3}
}  
```
---

## 📁 Project structure  
```
├── scripts/        # Implementation of algorithms
├── figs/           # Figures from repository
├── data/           # Example analysys results (small samples)
├── LICENSE
└── README.md
```
# 🚀 Methods

## 🌈 REACTIV (Rapid and EAsy Change detection in radar TIme-series by Variation coefficient)

REACTIV is an efficient and fast algorithm for detecting changes in radar (SAR) time series, designed for easy interpretation of results.  
The algorithm transforms each pixel into the HSV color space:

- **Hue** 🎨 → Time of maximum signal  
- **Saturation** 💡 → Intensity & reliability of the change  
- **Value** 🌟 → Maximum radar backscatter intensity  

After HSV → RGB conversion, the resulting visualization shows:  
- **Color** → Change timing  
- **Saturation** → Confidence level  
- **Brightness** → Strength of radar signal  

👉 Thresholding is applied to eliminate irrelevant fluctuations using empirically set parameters (typically `0.3` and `0.02`).

```js
// -------------------------------------------------------------
// Parameters: DATES, ASCENDING OR DESCENDING
var str2 = '2025-06-30';   // end date
var str1 = '2022-02-24';   // start date
//var str = 'DESCENDING';
var str = 'ASCENDING';

// -------------------------------------------------------------
// Thresholding parameters
var enableThresholding = true;   // Enable/disable thresholding
var cvThreshold = 0.3;           // Coefficient of variation threshold (0.75–1.0)
var intensityThreshold = -15;    // Maximum intensity threshold (dB)
var combinedThreshold = true;    // Use combined thresholds (true/false)
var croppalet = 0.6;             // Use <1 if you do not want the full HSV palette
```

📚 **Source**: [Reference Paper](https://www.mdpi.com/2072-4292/12/13/2089)  
💻 **Script**: [scripts/reactiv.js](./scripts/reactiv.js)

---

## 🌀 REACTIV Polarimetric

The polarimetric variant extends REACTIV by analyzing **VV** and **VH** polarization channels simultaneously.  

- Data is converted from logarithmic (dB) → linear amplitude.  
- Mean, standard deviation, and coefficient of variation are computed.  
- A covariance matrix is constructed, eigenvalues calculated, and the largest eigenvalue is normalized to capture polarimetric variability.  
- HSV encoding → RGB visualization reveals correlated changes across polarizations.

📚 **Source**: [Reference Paper](https://link.springer.com/article/10.1007/s12524-024-02005-x)  
💻 **Script**: [scripts/reactiv_polarimetric.js](./scripts/reactiv_polarimetric.js)

---

## ❄️ REACTIV FBG (Frozen Background Generation)

This variant estimates a **frozen background** (stable backscatter reference) by analyzing cumulative subsets of amplitude values and computing the coefficient of variation.  
It allows clearer distinction between **persistent structures** vs. **temporal anomalies**.

📚 **Source**: [Reference Paper](https://www.mdpi.com/2072-4292/12/11/1720)  
💻 **Script**: [scripts/reactiv_fbg.js](./scripts/reactiv_fbg.js)

---

## 🔍 Omnibus

The Omnibus algorithm applies a **likelihood ratio test** for covariance matrix homogeneity under the assumption of a Wishart distribution.  
- Detects subtle and gradual changes better than pairwise comparison.  
- Factorization enables pinpointing the exact change moment.

📚 **Source**: [Reference Paper](https://ieeexplore.ieee.org/document/7729878)  
💻 **Script**: [scripts/omnibus.js](./scripts/omnibus.js)

---

## 🏚️ PWTT (Pixel-Wise T-Test)

Developed for **conflict-related damage detection**:  
- Uses multi-temporal pixel standard deviation to reduce false alarms.  
- Requires a **pre-war reference period** + **short inference period**.  
- Produces probability maps of destruction in built-up areas.  

📚 **Source**: [Reference Paper](https://arxiv.org/abs/2405.06323)  
💻 **Script**: [scripts/pwtt.js](./scripts/pwtt.js)

---

## 🍂 Seasonality in REACTIV and Omnibus

<a href="./figs/seasonality.png" target="_blank">
  <img width="100%" src="./figs/seasonality.png" alt="Seasonality visualization">
</a>

To reduce false alarms from natural cycles, both algorithms are applied to **seasonally filtered data**:  
- REACTIV computes CV, intensity, and change timing only on seasonal subsets.

💻 **Script**: [scripts/reactiv_seasonal.js](./scripts/reactiv_seasonal.js)

- Omnibus compares covariance matrices across **same-season images**.

💻 **Script**: [scripts/omnibus_seasonal.js](./scripts/omnibus_seasonal.js)

---

# 📜 License

This repository is released under the **MIT License**.  
You are free to use, modify, and distribute the code for research and educational purposes, provided that proper attribution is given to the authors.  

See the [LICENSE](./LICENSE) file for full details.

---

# 📩 Contact

If you have questions about this repository, please contact the authors:

- **Kinga KARWOWSKA**, Military University of Technology, Department of Imagery Intelligence  
📧 email: [kinga.karwowska@wat.edu.pl](mailto:kinga.karwowska@wat.edu.pl)
- **Jakub SLESINSKI**, Military University of Technology, Department of Imagery Intelligence  
📧 email: [jakub.slesinski@wat.edu.pl](mailto:jakub.slesinski@wat.edu.pl)  

For issues and contributions, please use the [GitHub Issues](../../issues) section of this repository.  
We welcome collaborations and feedback from the research community.

