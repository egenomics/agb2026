#!/bin/bash -ue
python << 'PYTHON'
from PIL import Image, ImageDraw, ImageFont
from datetime import datetime

width = 300
height = 200

img = Image.new("RGB", (width, height), "white")
draw = ImageDraw.Draw(img)

try:
    font_title = ImageFont.truetype("DejaVuSansMono.ttf", 32)
    font_text  = ImageFont.truetype("DejaVuSansMono.ttf", 24)
except:
    font_title = ImageFont.load_default()
    font_text  = ImageFont.load_default()

lines = [
    "RAREFACTION SKIPPED",
    "",
    f"Sample count    : 58",
    f"Minimum required: 100",
    "",
    "Automatic rarefaction threshold selection was not",
    f"performed because the number of samples (58)",
    f"did not reach the minimum required (100).",
    "",
    "Action required: accumulate more samples and",
    "re-run the pipeline once the threshold is met."
]

y = 15

draw.text((20, y), lines[0], fill="black", font=font_title)
y += 20

for line in lines[2:]:
    draw.text((20, y), line, fill="black", font=font_text)
    y += 15

img.save("rarefaction_skipped_report.png")
PYTHON
