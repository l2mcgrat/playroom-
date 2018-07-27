import base64
import io
from PIL import Image
from sys import argv

picture = argv[1]

def pictureretriever(picture):

    if request.method == 'GET':
        if request.GET['pos'] and request.GET['frame']:
            pos = json.loads(request.GET['pos'])
            str_frame = json.loads(request.GET['frame'])

            # Converts the base64 string into a byte string, we need to encode
            # str_frame as utf-8 first otherwise python3.2 complains about unicode
            b64decoded_frame = base64.b64decode(str_frame.encode('utf-8'))

            # This puts the decoded image into a buffer so that I don't need to save
            # it to disk to use it in PIL
            byte_stream = io.BytesIO(b64decoded_frame)

            # Open the image and convert it to 8-bit greyscale (mono8)
            img = Image.open(byte_stream).convert('L')

            # Get the 8-bit pixel value
            pixel_val = img.getpixel((pos['x'], pos['y']))

            # Convert the 8-bit pixel value to 16-bit by holding the rations
            # i.e. n8 / (2^8 - 1) == x16 / (2^16 - 1)
            pixel_val = int(pixel_val / 255 * 65535)

    print pixel_val

    return

pictureretriever()
