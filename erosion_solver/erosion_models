#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import cgi, cgitb
import erosion_result

default_params = {"geometry": "jet",
                  "models": "finnie",
                  "vp": "1.0",
                  "alfa": "45.0",
                  "dp": "1e-6",
                  "rot": "8000.0",
                  "mp": "1.0",
                  "Hv": "200.0",
                  "rop": "3000.0",
                  "D": "0.1",
                  "mum": "0.001",
                  "rom": "1000.0",
                  "R": "5.0",
                  "D1": "0.1",
                  "D2": "0.05",
                  "BH": "600.0",
                  "P": "2.7e9",
                  "GF": "2.0",
                  "Fs": "0.5"}

def getParamsCGI():
    form = cgi.FieldStorage()
    if not form.keys():
        return default_params
    params={}
    params["geometry"] = form.getfirst("geometry", "")
    params["models"] = form.getfirst("models","")
    params["vp"] = form.getfirst("vp", "")
    params["alfa"] = form.getfirst("alfa", "")
    params["dp"] = form.getfirst("dp", "")
    params["rot"] = form.getfirst("rot", "")
    params["mp"] = form.getfirst("mp", "")
    params["Hv"] = form.getfirst("Hv", "")
    params["rop"] = form.getfirst("rop", "")
    params["D"] = form.getfirst("D", "")
    params["mum"] = form.getfirst("mum", "")
    params["rom"] = form.getfirst("rom", "")
    params["R"] = form.getfirst("R", "")
    params["D1"] = form.getfirst("D1", "")
    params["D2"] = form.getfirst("D2", "")
    params["BH"] = form.getfirst("BH", "")
    params["P"] = form.getfirst("P", "")
    params["GF"] = form.getfirst("GF", "")
    params["Fs"] = form.getfirst("Fs", "")

    return params

def main():
    print "Content-Type: text/html\n\n"
    paramst = getParamsCGI()
    model = loadHTMLModelHeaderFooter("erosion_models.html","static//header_lemt","static//footer_lemt")
    convertErrorMessage = "Parâmetro inválido"
    modelErrorMessage = "Modelo inválido"
    readError = False
    results = [""] * 6
    list_select = {"finnie": ["vp", "alfa", "rot", "mp", "P"],
                   "oka": ["vp", "alfa", "dp", "rot", "mp", "Hv"],
                   "ahlert": ["vp", "alfa", "rot", "mp", "BH", "Fs"],
                   "dnv": ["vp", "alfa", "rot", "mp"],
                   "haugen": ["vp", "alfa", "rot", "mp"],
                   "huang": ["vp", "alfa", "dp", "rot", "mp", "rop"],
                   "arabnejad": ["vp", "alfa","rot", "mp", "Fs"],
                   "pipe": ["vp", "mp", "D"],
                   "joint": ["vp", "alfa", "dp", "rot", "mp", "D", "rom"],
                   "bend": ["vp", "dp", "rot", "mp", "rop", "D", "mum", "rom", "R", "GF"],
                   "tee": ["vp", "dp", "rot", "mp", "rop", "D", "mum", "rom", "GF"],
                   "reducer": ["vp", "alfa", "dp", "rot", "mp","rom", "D1", "D2", "GF"]}
    geometry_used = paramst["geometry"]
    model_used = paramst["models"]
    params_used = []
    if geometry_used == "jet":
        index_params_used = list_select[model_used]
    else:
        index_params_used = list_select[geometry_used]

    for i in index_params_used:
        try:
            params_used.append(float(paramst[i]))
        except:
            paramst[i] = convertErrorMessage
            readError = True
    if not readError:
        try:
            results = erosion_result.models(geometry_used, model_used, params_used)
            results = ["{0:.4g}".format(x) for x in results]
        except:
            results = [modelErrorMessage] * 2
    print generateHTML(model,paramst,results)
    return 0

def loadHTMLModelHeaderFooter(filename, headerfile, footerfile):
    headerplaceholder = "##header##"
    footerplaceholder = "##footer##"

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), filename)
    with open(path, 'r') as f:
        model = f.read()
    path_header = os.path.join(os.path.dirname(os.path.realpath(__file__)), headerfile)
    if os.path.isfile(path_header):
	with open(path_header, 'r') as f:
            header = f.read()
    else:
        header = ""
    path_footer = os.path.join(os.path.dirname(os.path.realpath(__file__)), footerfile)
    if os.path.isfile(path_footer):
        with open(path_footer, 'r') as f:
            footer = f.read()
    else:
        footer = ""

    model = model.replace(headerplaceholder,header)
    model = model.replace(footerplaceholder,footer)
    return model

def generateHTML(model,params,results):
    geo = params["geometry"]
    mod = params["models"]
    new_params = params
    newHTML = model
    patternPar = "##{0}p##"
    patternRes = "##{0}r##"
    patternMod = "##{0}m##"
    patternGeo = "##{0}g##"
    for i, value in params.items():
        newHTML = newHTML.replace(patternPar.format(i),value)
    for i in range(len(results)):
        newHTML = newHTML.replace(patternRes.format(str(i)),results[i])
    newHTML = newHTML.replace(patternGeo.format(str(geo)), "selected")
    newHTML = newHTML.replace(patternMod.format(str(mod)), "selected")

    return newHTML

if __name__ == "__main__":
    sys.exit(main())
