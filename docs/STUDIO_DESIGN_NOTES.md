# Studio 1.0.0 design notes

DockRIFT Studio deliberately follows a product-style scientific-workbench model rather than a generic form-based dashboard. The interface is custom HTML/CSS/JavaScript served locally by Python. Plotly JavaScript is served from the installed package, so the core visualization layer does not require a CDN.

The visual system uses a deep navy scientific workspace, electric cyan/blue/violet accents, glass-like panels, a persistent navigation rail, responsive metric cards, censor-aware badges, and information-dense but non-cluttered charts.

The software architecture preserves one crucial rule: the GUI is presentation/orchestration only. Numerical definitions remain in `ranks.py`, `metrics.py`, `rrl.py`, `bootstrap.py`, and `validation.py`. The browser does not contain a second implementation of DockRIFT equations.
