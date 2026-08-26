"""Compatibility launcher for DockRIFT Studio."""
from .gui.server import main, serve
__all__=["main","serve"]
if __name__ == "__main__": main()
