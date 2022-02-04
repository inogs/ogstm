#! /usr/bin/python3
from xml.dom import minidom
import sys

def file2stringlist(filename):
    LIST=[]
    filein=open(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST
xmldoc = minidom.parse('BFMtab.xml')

ORIG_NAMELIST=file2stringlist("namelist.passivetrc")
endnamelist1 = ORIG_NAMELIST.index("/")
endnamelist2 = ORIG_NAMELIST.index("/",endnamelist1+1)
endnamelist3 = ORIG_NAMELIST.index("/",endnamelist2+1)


XML_MODELVARS={}
XML_DIA__VARS={}
XML_DIA__DUMP={}
XML_DIA2dVARS={}
XML_DIA2dDUMP={}

NODES=xmldoc.getElementsByTagName("ModelVars")[0].getElementsByTagName("var")
for node in NODES:
    var = str(node.getAttribute("name"))
    max = float(node.getAttribute("maxvalue"))
    high= node.getAttribute("highfreq")=='true'
    if var not in XML_MODELVARS.keys(): # check about BFMtab formatting
        XML_MODELVARS[var] = [max,high]
    else:
        print(var + " is already defined.")
        sys.exit(1)

NODES=xmldoc.getElementsByTagName("Diagnostics")[0].getElementsByTagName("var")
for node in NODES:
    var = str(node.getAttribute("name"))
    high= node.getAttribute("highfreq")=='true'
    dump= node.getAttribute("dump")=='true'
    if var not in XML_DIA__VARS.keys():
        XML_DIA__VARS[var] = high
    else:
        print(var + " is already defined.")
        sys.exit(1)
    if var not in XML_DIA__DUMP.keys():
        XML_DIA__DUMP[var] = dump
    else:
        print(var + " is already defined.")
        sys.exit(1)  

NODES=xmldoc.getElementsByTagName("Diagnostics_2D")[0].getElementsByTagName("var")
for node in NODES:
    var = str(node.getAttribute("name"))
    high= node.getAttribute("highfreq")=='true'
    dump= node.getAttribute("dump")=='true'
    if var not in XML_DIA2dVARS.keys():
        XML_DIA2dVARS[var] = high
    else:
        print(var + " is already defined.")
        sys.exit(1)
    if var not in XML_DIA2dDUMP.keys():
        XML_DIA2dDUMP[var] = dump
    else:
        print(var + " is already defined.")
        sys.exit(1) 
        



MODEL_HF=[]
MODEL_LF=[]
DIA_HF  =[]
DIA_LF  =[]

NAMELIST_NEW=[]
NAMELIST_NEW.append("&NATTRC \n")

SECTION_NAMELIST=ORIG_NAMELIST[:endnamelist1]
for il, line in enumerate(SECTION_NAMELIST):
    if line.find("ctrcnm") != -1:
        quote_1=line.find("\"")
        quote_2=line.find("\"",quote_1+1)
        varname=line[quote_1+1:quote_2]
        if varname in XML_MODELVARS.keys():
            par_1 = line.find("(")
            par_2 = line.find(")")
            ind   = int(line[par_1+1:par_2])
            
            NAMELIST_NEW.append(line +"\n")
            NAMELIST_NEW.append(SECTION_NAMELIST[il+1]+"\n")
            NAMELIST_NEW.append("   ctrmax(%d)=%e \n" % (ind, XML_MODELVARS[varname][0]))
            NAMELIST_NEW.append("   ctr_hf(%d)=%d\n" % (ind, XML_MODELVARS[varname][1]))
            NAMELIST_NEW.append("\n")
            
            if XML_MODELVARS[varname][1]:
                MODEL_HF.append(varname)
            MODEL_LF.append(varname)
        else:
            print("Error: " + varname + " not defined !")
            sys.exit(1)
            
NAMELIST_NEW.append("/\n\n")

SECTION_NAMELIST=ORIG_NAMELIST[endnamelist1:endnamelist2]
NAMELIST_NEW.append("&NATTRC_DIAG\n")
ivar=1
for il, line in enumerate(SECTION_NAMELIST):
    if line.find("dianm") != -1:
        quote_1=line.find("\"")
        quote_2=line.find("\"",quote_1+1)
        varname=line[quote_1+1:quote_2]
        if varname in XML_DIA__VARS.keys():
            par_1 = line.find("(")
            par_2 = line.find(")")

            _,right_side=line.rsplit("=")
            NAMELIST_NEW.append("   dianm(%d)=%s\n" %(ivar,right_side))
            _,right_side = SECTION_NAMELIST[il+1].rsplit("=")
            NAMELIST_NEW.append("   diaun(%d)=%s\n" %(ivar,right_side))
            NAMELIST_NEW.append("   diahf(%d)=%d\n" % (ivar, XML_DIA__VARS[varname]))
            NAMELIST_NEW.append("   diaWR(%d)=%d\n" % (ivar, XML_DIA__DUMP[varname]))

            NAMELIST_NEW.append("\n")
            ivar=ivar+1

            if XML_DIA__VARS[varname]:
                if XML_DIA__DUMP[varname]: DIA_HF.append(varname)
            if XML_DIA__DUMP[varname]: DIA_LF.append(varname)
        else:
            print("Error: " + varname + " not defined !")
            sys.exit(1)
NAMELIST_NEW.append("/ \n\n")

SECTION_NAMELIST=ORIG_NAMELIST[endnamelist2:]
NAMELIST_NEW.append("&NATTRC_DIAG_2D\n")
ivar=1
for il, line in enumerate(SECTION_NAMELIST):    
    if line.find("dianm") != -1:
        quote_1=line.find("\"")
        quote_2=line.find("\"",quote_1+1)
        varname=line[quote_1+1:quote_2]
        if varname in XML_DIA2dVARS.keys():
            par_1 = line.find("(")
            par_2 = line.find(")")

            _,right_side=line.rsplit("=")
            NAMELIST_NEW.append("   dianm_2d(%d)=%s\n" %(ivar,right_side))
            _,right_side = SECTION_NAMELIST[il+1].rsplit("=")
            NAMELIST_NEW.append("   diaun_2d(%d)=%s\n" %(ivar,right_side))

            NAMELIST_NEW.append("   diahf_2d(%d)=%d\n" % (ivar, XML_DIA2dVARS[varname]))
            NAMELIST_NEW.append("   diaWR_2d(%d)=%d\n" % (ivar, XML_DIA2dDUMP[varname]))

            NAMELIST_NEW.append("\n")
            ivar=ivar+1

            if XML_DIA2dVARS[varname]:
                if XML_DIA2dDUMP[varname] : DIA_HF.append(varname)
            if XML_DIA2dDUMP[varname] :DIA_LF.append(varname)
        else:
            print("Error: " + varname + " not defined !")
            sys.exit(1)
NAMELIST_NEW.append("/\n\n")



      

F=open("namelist.passivetrc_new","w")
F.writelines(NAMELIST_NEW)
F.close()

DUMP_DIA=[var for var in XML_DIA__DUMP.keys() if XML_DIA__DUMP[var] ]



printout=True
if printout:
    print("HF variables : " , len(MODEL_HF), "State", len(DIA_HF), "Diagnostics")
    print("LF variables : " , len(MODEL_LF), "State", len(DIA_LF), "Diagnostics", len(DUMP_DIA), "Dumped")

    print("HF LIST:")
    for var in MODEL_HF: print(var)
    for var in DIA_HF : print(var)

    print("\n\nDIA DUMPED:")
    for var in DUMP_DIA : print("  ", var)
