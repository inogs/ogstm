from xml.dom import minidom
import sys

def file2stringlist(filename):
    LIST=[]
    filein=file(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST
xmldoc = minidom.parse('BFMtab.ben.xml')

ORIG_NAMELIST=file2stringlist("namelist.ben")
endnamelist1 = ORIG_NAMELIST.index("/")
endnamelist2 = ORIG_NAMELIST.index("/",endnamelist1+1)
endnamelist3 = ORIG_NAMELIST.index("/",endnamelist2+1)

XML_MODELVARS={}
XML_DIA__VARS={}
XML_DIA2dVARS={}

NODES=xmldoc.getElementsByTagName("ModelVars")[0].getElementsByTagName("var")
for node in NODES:
    var = str(node.getAttribute("name"))
    max = float(node.getAttribute("maxvalue"))
    high= node.getAttribute("highfreq")=='true'
    if not XML_MODELVARS.has_key(var): # check about BFMtab formatting
        XML_MODELVARS[var] = [max,high]
    else:
        print var + " is already defined."
        sys.exit(1)

NODES=xmldoc.getElementsByTagName("Diagnostics")[0].getElementsByTagName("var")
for node in NODES:
    var = str(node.getAttribute("name"))
    high= node.getAttribute("highfreq")=='true'
    if not XML_DIA__VARS.has_key(var):
        XML_DIA__VARS[var] = high
    else:
        print var + " is already defined."
        sys.exit(1)        

NODES=xmldoc.getElementsByTagName("Diagnostics_2D")[0].getElementsByTagName("var")
for node in NODES:
    var = str(node.getAttribute("name"))
    high= node.getAttribute("highfreq")=='true'
    if not XML_DIA2dVARS.has_key(var):
        XML_DIA2dVARS[var] = high
    else:
        print var + " is already defined."
        sys.exit(1) 



MODEL_HF=[]
MODEL_LF=[]
DIA_HF  =[]
DIA_LF  =[]

NAMELIST_NEW=[]
NAMELIST_NEW.append("&NATBEN \n")

SECTION_NAMELIST=ORIG_NAMELIST[:endnamelist1]
for il, line in enumerate(SECTION_NAMELIST):
    if line.find("cbennm") != -1:
        quote_1=line.find("\"")
        quote_2=line.find("\"",quote_1+1)
        varname=line[quote_1+1:quote_2]
        if XML_MODELVARS.has_key(varname):
            par_1 = line.find("(")
            par_2 = line.find(")")
            ind   = int(line[par_1+1:par_2])
            
            NAMELIST_NEW.append(line +"\n")
            NAMELIST_NEW.append(SECTION_NAMELIST[il+1]+"\n")
            NAMELIST_NEW.append("   cbemax(%d)=%e \n" % (ind, XML_MODELVARS[varname][0]))
            NAMELIST_NEW.append("   cbe_hf(%d)=%d \n" % (ind, XML_MODELVARS[varname][1]))
            NAMELIST_NEW.append("\n")
            
            if XML_MODELVARS[varname][1]:
                MODEL_HF.append(varname)
            else:                                
                MODEL_LF.append(varname)
        else:
            print "Error: " + varname + " not defined !"
            sys.exit(1)
            
NAMELIST_NEW.append("/ \n\n")

SECTION_NAMELIST=ORIG_NAMELIST[endnamelist1:endnamelist2]
NAMELIST_NEW.append("&NATBEN_DIAG\n")
for il, line in enumerate(SECTION_NAMELIST):
    if line.find("dibnm") != -1:
        quote_1=line.find("\"")
        quote_2=line.find("\"",quote_1+1)
        varname=line[quote_1+1:quote_2]
        if XML_DIA__VARS.has_key(varname):
            par_1 = line.find("(")
            par_2 = line.find(")")
            ind   = int(line[par_1+1:par_2])
            
            NAMELIST_NEW.append(line +"\n")
            NAMELIST_NEW.append(SECTION_NAMELIST[il+1]+"\n")
            NAMELIST_NEW.append("   dibhf(%d)=%d \n" % (ind, XML_DIA__VARS[varname]))
            NAMELIST_NEW.append("\n")
            
            if XML_DIA__VARS[varname]:
                DIA_HF.append(varname)
            else:
                DIA_LF.append(varname)
        else:
            print "Error: " + varname + " not defined !"
            sys.exit(1)
NAMELIST_NEW.append("/ \n\n")

SECTION_NAMELIST=ORIG_NAMELIST[endnamelist2:]
NAMELIST_NEW.append("&NATBEN_DIAG_2D\n")
for il, line in enumerate(SECTION_NAMELIST):    
    if line.find("dibnm") != -1:
        quote_1=line.find("\"")
        quote_2=line.find("\"",quote_1+1)
        varname=line[quote_1+1:quote_2]
        if XML_DIA2dVARS.has_key(varname):
            par_1 = line.find("(")
            par_2 = line.find(")")
            ind   = int(line[par_1+1:par_2])
            
            NAMELIST_NEW.append(line.replace("dibnm","dibnm_2d")+"\n" )
            NAMELIST_NEW.append(SECTION_NAMELIST[il+1].replace("dibun","dibun_2d") +"\n")
            NAMELIST_NEW.append("   dibhf_2d(%d)=%d \n" % (ind, XML_DIA2dVARS[varname]))
            NAMELIST_NEW.append("\n")
            
            if XML_DIA2dVARS[varname]:
                DIA_HF.append(varname)
            else:
                DIA_LF.append(varname)
        else:
            print "Error: " + varname + " not defined !"
            sys.exit(1)
NAMELIST_NEW.append("/ \n\n")



      

F=open("namelist.ben_new","w")
F.writelines(NAMELIST_NEW)
F.close()

printout=True
if printout:
    print "HF variables : " , len(MODEL_HF), "Model", len(DIA_HF), "Diagnostics"
    print "LF variables : " , len(MODEL_LF), "Model", len(DIA_LF), "Diagnostics"
    
    print "HF LIST:"
    for var in MODEL_HF: print var
    for var in DIA_HF : print var
