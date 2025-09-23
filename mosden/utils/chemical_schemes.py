def MSBR_scheme(include_long=True,
                rate_scaling=1) -> dict[str: float]:
    rate = 1/20 * rate_scaling
    repr_dict = {'Kr': rate,
                'Xe': rate,
                'Se': rate,
                'Nb': rate,
                'Mo': rate,
                'Tc': rate,
                'Ru': rate,
                'Rh': rate,
                'Pd': rate,
                'Ag': rate,
                'Sb': rate,
                'Te': rate,
                } 
    if not include_long:
        return repr_dict

    rate = 1/(3*24*3600) * rate_scaling
    more_data = {'Pa': rate}
    repr_dict.update(more_data)

    rate = 1/(50*24*3600) * rate_scaling
    more_data = {'Y': rate,
                'La': rate,
                'Ce': rate,
                'Pr': rate,
                'Nd': rate,
                'Pm': rate,
                'Sm': rate,
                'Gd': rate,
    }
    repr_dict.update(more_data)

    rate = 1/(60*24*3600) * rate_scaling
    more_data = {'Br': rate,
                'I': rate
    }
    repr_dict.update(more_data)

    rate = 1/(200*24*3600) * rate_scaling
    more_data = {'Zr': rate,
                'Cd': rate,
                'In': rate,
                'Sn': rate
    }
    repr_dict.update(more_data)

    rate = 1/(500*24*3600) * rate_scaling
    more_data = {'Eu': rate}
    repr_dict.update(more_data)

    return repr_dict