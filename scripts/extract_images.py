import os
import json
import base64
import sys

def extract_images_from_notebook(notebook_path, output_dir="./"):
    # Open the notebook file
    with open(notebook_path, 'r') as f:
        notebook_data = json.load(f)

    image_count = 0

    # Loop through each cell in the notebook
    for cell in notebook_data.get("cells", []):
        # Check if cell has outputs
        if "outputs" in cell:
            for output in cell["outputs"]:
                # Check if output contains mage data
                if "data" in output and "image/png" in output["data"]:
                    # Decode the base64 image data
                    image_data = output["data"]["image/png"]
                    image_bytes = base64.b64decode(image_data)

                    # Save the image as a PNG file
                    image_path = os.path.join(output_dir, "image_{}.png".format(image_count))
                    with open(image_path, 'wb') as img_file:
                        img_file.write(image_bytes)

                    image_count += 1

    print("Extracted {} images to {}".format(image_count,output_dir))

# Use the function
file = sys.argv[1]
extract_images_from_notebook(file)
