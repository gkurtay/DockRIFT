# DockRIFT 1.0.0 architecture

The numerical scientific core is implemented in `dockrift.rrl`, `dockrift.ranks`, `dockrift.metrics`, and `dockrift.validation`. `dockrift.gui.logic` is an adapter that calls those functions and prepares table/plot/report payloads. `dockrift.gui.server` is a thin local HTTP API and static-file server. `dockrift/gui/static/` contains the custom browser presentation layer.

This separation is intentional: GUI presentation does not define an alternative numerical method. CLI, Studio, and regression tests use the same core functions.

## Local processing flow

score table → normalization/QC → in-process pandas DataFrame → core metrics → JSON/Plotly payload → browser visualization → local reproducibility export

DockRIFT requires no cloud database, telemetry endpoint, or external analytics service.
