from django import template
import math

register = template.Library()


@register.filter
def sci(value):
    """
    Format a float in scientific notation like  2.3×10⁻³⁹
    Usage: {{ result.evalue|sci }}
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return value

    if value == 0:
        return "0"

    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10 ** exponent)

    # Unicode superscript digits and minus sign
    sup_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    sup_exp = str(exponent).translate(sup_map)

    return f"{mantissa:.1f}×10{sup_exp}"
