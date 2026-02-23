# Flask has been removed; reports are now rendered via Jinja2 directly.
# This script is kept for reference but no longer starts a Flask server.
if __name__ == "__main__":
    raise SystemExit(
        "The Flask server has been replaced by direct Jinja2 rendering. "
        "Use the Reporter class to generate HTML reports instead."
    )
