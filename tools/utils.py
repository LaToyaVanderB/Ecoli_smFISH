import numpy as np
from copy import copy

# with eternal thanks to the interwebs:
def translate_image(img, tx, ty):
    N, M = img.shape
    img_translated = np.zeros_like(img)

    img_translated[max(tx, 0):M + min(tx, 0), max(ty, 0):N + min(ty, 0)] = img[-min(tx, 0):M - max(tx, 0),
                                                                             -min(ty, 0):N - max(ty, 0)]

    img_translated_cropped = img[-min(tx, 0):M - max(tx, 0), -min(ty, 0):N - max(ty, 0)]
    return img_translated, img_translated_cropped

