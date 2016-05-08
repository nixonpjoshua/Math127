from common_lines import common_lines
from vanheel import estimate_orientations


def reconstruct_orientations(images, num_lines, granularity):
    """
    :param images: the images which we want to find orientations for
    :return: the estimated orientation of every images
    """
    return estimate_orientations(common_lines(images, num_lines, granularity))
