# Literature LLE Data for Binary Solvent Systems

Compiled from IUPAC-NIST Solubility Database and peer-reviewed publications for comparison with openCOSMO-RS calculations.

---

## 1. Water + 1-Butanol

### Mutual Solubility at 298.2 K (25°C) - IUPAC Recommended

| Phase | Mole Fraction 1-Butanol | Mole Fraction Water |
|-------|-------------------------|---------------------|
| Water-rich | 0.0191 | 0.9809 |
| Butanol-rich | 0.488 | 0.512 |

### Critical Solution Temperature
- **UCST**: ~120-125°C (393-398 K)
- **Critical composition**: x_butanol ≈ 0.11

### Key References
- IUPAC Solubility Data Series, Volume 15 - "Alcohols with Water"
- IUPAC-NIST Solubility Data Series 82 - Alcohols with Water (Revised)
- Hill and Malisoff (1926), J. Am. Chem. Soc. 48, 918-927
- Moriyoshi et al. (1978), J. Chem. Thermodyn. 10, 1003

### URLs
- https://srdata.nist.gov/solubility/
- https://pubs.acs.org/doi/10.1021/ja01415a011

---

## 2. Water + n-Hexane

### Mutual Solubility Data

**At 298 K (25°C):**
| Phase | Mole Fraction | Concentration |
|-------|---------------|---------------|
| Water in hexane-rich | ~4.8 × 10⁻⁴ | ~0.0036 M |
| Hexane in water-rich | ~2.3 × 10⁻⁶ | ~9.8 mg/L |

**Temperature Correlation (293-353 K):**
```
s (10² g water/100g sln) = 86.5345 - 0.6183*T + 0.00113*T²
```

### Thermodynamic Parameters
- Enthalpy of solution (water in hexane): ΔH_sln = 28.7 kJ/mol
- Heat capacity of solution: ΔCp_sln = 73 J/(K·mol)

### Key References
- IUPAC-NIST Solubility Data Series, Volume 37
- Maczynski et al. (2005), J. Phys. Chem. Ref. Data, 34(2), 709-753
- Polak & Lu (1973), Can. J. Chem., 51, 4018-4023
- Tsonopoulos & Wilson (1983), AIChE Journal, 29, 990-999
- Dupeux et al. (2019), Fluid Phase Equilibria

### URLs
- https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=37_282
- https://pubs.aip.org/aip/jpr/article-abstract/34/2/709/242032

---

## 3. Water + Ethyl Acetate

### Mutual Solubility at 20°C (293.15 K)

| Phase | Composition |
|-------|-------------|
| Water-rich | 8.7 wt% ethyl acetate |
| Ethyl acetate-rich | 3.3 wt% water |

### At 25°C (298.15 K)
- Ethyl acetate in water: ~8.0-8.3 wt% (~83-87 g/L)
- Water in ethyl acetate: ~3.0 wt%

### Heterogeneous Azeotrope
- Boiling point: 70.4°C
- Composition: ~69-70 wt% ethyl acetate

### Key References
- IUPAC-NIST Solubility Data Series, Volume 88
- Trofimova et al. (2020), Fluid Phase Equilibria
- Toikka et al. (2012), Fluid Phase Equilibria

### URLs
- https://pubs.aip.org/aip/jpr/article/38/4/1093/1059855
- https://www.sciencedirect.com/science/article/abs/pii/S0378381219303826
- https://macro.lsu.edu/howto/solvents/ethylacetate.htm

---

## 4. Water + Tetrahydrofuran (THF)

### Unusual Closed-Loop Miscibility Gap (Type VI)

| Property | Value |
|----------|-------|
| LCST | ~71.8°C |
| UCST | ~137.1°C |
| Immiscibility range | 0.28-0.72 weight fraction THF |

### Azeotrope Data
- Temperature: 63.4-64°C at 101.325 kPa
- Composition: 0.82 mole fraction THF (~93.3 wt%)

### Key References
- Matous et al. (1972), Collect. Czech. Chem. Commun. 37, 2653-2663
- Kiyohara & Benson (1977), Can. J. Chem. 55(8), 1354-1359
- Hayduk et al. (1973), J. Chem. Eng. Data 18, 373-376

### URLs
- https://cdnsciencepub.com/doi/10.1139/v77-187
- https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5b09770
- https://webbook.nist.gov/cgi/cbook.cgi?ID=C109999

---

## 5. Methanol + n-Hexane

### Upper Critical Solution Temperature
| Source | UCST |
|--------|------|
| IUPAC-NIST Recommended | 306.8 K (33.7°C) |
| Hradetzky & Lempe (1991) | ~308 K (~35°C) |
| Trejo et al. (2006) | 308.3 K (35.2°C) |

### Binodal Curve Data at 293.2 K (20°C)

| Phase | Methanol Mole Fraction | Approx. wt% |
|-------|------------------------|-------------|
| Hexane-rich | 0.210 | ~8.4% |
| Methanol-rich | 0.822 | ~62.5% |

### Temperature Dependence (approximate)

| T (K) | x_MeOH (hexane-rich) | x_MeOH (methanol-rich) |
|-------|----------------------|------------------------|
| 273 | ~0.15 | ~0.87 |
| 283 | ~0.18 | ~0.85 |
| 293 | 0.21 | 0.82 |
| 298 | ~0.24 | ~0.78 |
| 303 | ~0.30 | ~0.70 |
| 305.7 | ~0.35 | ~0.35 (near UCST) |

### Key References
- Trejo et al. (2006), J. Chem. Eng. Data 51, 1070-1075
- Hradetzky & Lempe (1991), Fluid Phase Equilibria 69, 285-301
- McLure et al. (1997), J. Chem. Soc. Faraday Trans.
- IUPAC-NIST Solubility Data Series Vol. 56

### URLs
- https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=69_21
- https://pubs.acs.org/doi/10.1021/je0505321
- https://www.sciencedirect.com/science/article/abs/pii/037838129190040E

---

## 6. Acetonitrile + n-Hexane

### Critical Solution Temperature
- **UCST**: 350.2 K (77.05°C)
- **Critical mole fraction acetonitrile**: x_c = 0.583

### Binodal Curve Data

| T (K) | x_ACN (hexane-rich) | x_ACN (ACN-rich) |
|-------|---------------------|------------------|
| 283.2 | 0.035 | 0.970 |
| 293 | ~0.05 | ~0.95 |
| 298 | ~0.06 | ~0.93 |
| 303 | ~0.07 | ~0.91 |
| 313 | ~0.10 | ~0.87 |
| 323 | ~0.14 | ~0.82 |
| 333 | ~0.21 | ~0.76 |
| 343 | ~0.32 | ~0.70 |
| 350.0 | 0.503 | 0.657 |
| 350.2 (UCST) | 0.583 | 0.583 |

### Key References
- Sinegubova (1978) - Primary data source
- McLure et al. (1982), Fluid Phase Equilibria 8, 271
- IUPAC-NIST Solubility Data Series 78
- J. Chem. Eng. Data (2024) - Recent update

### URLs
- https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=78_61
- https://pubs.acs.org/doi/10.1021/acs.jced.4c00104
- https://www.researchgate.net/publication/238956506

---

## Summary Table: Key Systems for Validation

| System | UCST/LCST | x₁ (phase 1) @ 25°C | x₁ (phase 2) @ 25°C |
|--------|-----------|---------------------|---------------------|
| Water + 1-Butanol | UCST ~398 K | 0.019 | 0.488 |
| Water + n-Hexane | No CST (immiscible) | 4.8×10⁻⁴ | 2.3×10⁻⁶ |
| Water + Ethyl Acetate | No CST | ~0.05 (wt: 8.7%) | ~0.02 (wt: 3.3%) |
| Water + THF | LCST 345K, UCST 410K | Miscible at 25°C | Miscible at 25°C |
| Methanol + n-Hexane | UCST ~307 K | 0.21 | 0.82 |
| Acetonitrile + n-Hexane | UCST 350.2 K | 0.06 | 0.93 |

---

## Primary Database Sources

1. **IUPAC-NIST Solubility Database**: https://srdata.nist.gov/solubility/
2. **Dortmund Data Bank (DDB)**: https://www.ddbst.com
3. **DECHEMA Chemistry Data Series**: https://dechema.de
4. **NIST Chemistry WebBook**: https://webbook.nist.gov

---

*Compiled: January 2026*
*For comparison with openCOSMO-RS calculations*
