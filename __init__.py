import _pywfl
import filters
import fit


def factorial(n):
    """ Get :math:`n!=n \dot (n-1) \cdots 3 \dot 2 \dot 1`

    Args:
        n (int): number to be calculated factorial.

    Returns:
        float : get :math:`n!=n \dot (n-1) \cdots 3 \dot 2 \dot 1`

    Notes:
        Use recursive method to do factorial.
        FACTORIAL is slower than WFL_FACTORIAL when call it many times.


"""
    return _pywfl.py_factorial(n)
