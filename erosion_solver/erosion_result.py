# -*- coding: utf-8 -*-
from erosion_equations import finnie, oka, ahlert, dnv, haugen, huang, arabnejad
from erosion_equations import smooth, welded, bend, blinded, reducer

def models(geometry, model, params):

    model_dict = {"finnie": finnie, "oka": oka, "ahlert": ahlert, "dnv": dnv, "haugen": haugen, "huang": huang, "arabnejad": arabnejad}
    geometry_dict = {"pipe": smooth, "joint": welded, "bend": bend, "tee": blinded, "reducer": reducer}
    if geometry == "jet":
        if model in model_dict:
            func = model_dict[model]
            return func(*params)
    else:
        if geometry in geometry_dict:
            func = geometry_dict[geometry]
            return func(*params)

    # TODO: error checking
