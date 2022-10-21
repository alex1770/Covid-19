import classifycog, classifygisaid
from variantaliases import aliases

ecache={}
def expandlin(lin):
  if lin in ecache: return ecache[lin]
  for (short,long) in aliases:
    s=len(short)
    if lin[:s+1]==short+".": ecache[lin]=long+lin[s:];return ecache[lin]
  ecache[lin]=lin
  return lin

ccache={}
def contractlin(lin):
  if lin in ccache: return ccache[lin]
  lin=expandlin(lin)
  for (short,long) in aliases:
    l=len(long)
    if lin[:l+1]==long+".": ccache[lin]=short+lin[l:];return ccache[lin]
  ccache[lin]=lin
  return lin

# Manual addition to automatic classification, because the official classifications are so far behind that the learner can't be trained properly
def manualclassifycog(cl,mutations):

  # Manually classify BA.2.75*
  if ("S:G446S" in mutations)+("S:G257S" in mutations)+("S:K147E" in mutations)+("S:G339H" in mutations)>=2: # BA.2.75*
    if "S:D574V" in mutations:
      if "S:R346T" in mutations: cl="BL.1"
      else: cl="BA.2.75.1"
    elif "S:R346T" in mutations and "S:F486S" in mutations and "S:D1199N" in mutations: cl="BA.2.75.2"
    elif "orf1ab:S1221L" in mutations and "orf1ab:V6107I" in mutations: cl="BA.2.75.3"
    elif "S:L452R" in mutations: cl="BA.2.75.4"
    elif "S:K356T" in mutations:
      if "S:F490S" in mutations and "S:R346T" in mutations: cl="BN.1"
      else: cl="BA.2.75.5"
    else: cl="BA.2.75"
  
  # BE.1.1.stuff -> BQ.stuff
  if cl=="BE.1.1" and "S:K444T" in mutations: cl="BE.1.1.1"
  if cl=="BE.1.1.1" and "S:N460K" in mutations: cl="BQ.1"
  if cl=="BQ.1":
    if "S:R346T" in mutations:
      if "orf1ab:V6040A" in mutations: cl="BQ.1.1.1"
      else: cl="BQ.1.1"
    elif "S:I666V" in mutations: cl="BQ.1.2"
    elif "S:E619Q" in mutations: cl="BQ.1.3"
    elif "S:R190T" in mutations: cl="BQ.1.4"
    elif "orf1ab:T5477I" in mutations: cl="BQ.1.5"
    elif "orf1ab:S505F" in mutations: cl="BQ.1.8"

  # Offshoots of BA.2.75.*
  if cl=="BA.2.75.2" and "S:T604I" in mutations and "S:L452R" in mutations: cl="CA.1"
  
  if cl=="BA.2.75.3" and "S:F486S" in mutations: cl="BM.1"
  if cl=="BM.1" and "S:R346T" in mutations: cl="BM.1.1"
  if cl=="BM.1.1" and "S:F490S" in mutations: cl="BM.1.1.1"

  if cl=="BA.2.75.4" and "S:R346T" in mutations and "S:F486I" in mutations: cl="BR.2"

  if cl=="BA.2.3" and "S:E484R" in mutations: cl="BA.2.3.20"
  
  # Recombinant
  if "S:G252V" in mutations and ("S:F486S" in mutations)+("S:F490S" in mutations)+("S:R346T" in mutations)>=2: cl="XBB.1"
  
  return cl

# Manual addition to automatic classification, because the official classifications are so far behind that the learner can't be trained properly
def manualclassifygisaid(cl,mutations):

  # Manually classify BA.2.75*
  if ("Spike_G446S" in mutations)+("Spike_G257S" in mutations)+("Spike_K147E" in mutations)+("Spike_G339H" in mutations)>=2: # BA.2.75*
    if "Spike_D574V" in mutations:
      if "Spike_R346T" in mutations: cl="BL.1"
      else: cl="BA.2.75.1"
    elif "Spike_R346T" in mutations and "Spike_F486S" in mutations and "Spike_D1199N" in mutations: cl="BA.2.75.2"
    elif "NSP3_S403L" in mutations and "NSP14_V182I" in mutations: cl="BA.2.75.3"
    elif "Spike_L452R" in mutations: cl="BA.2.75.4"
    elif "Spike_K356T" in mutations:
      if "Spike_F490S" in mutations and "Spike_R346T" in mutations: cl="BN.1"
      else: cl="BA.2.75.5"
    else: cl="BA.2.75"
  
  # BE.1.1.stuff -> BQ.stuff
  if cl=="BE.1.1" and "Spike_K444T" in mutations: cl="BE.1.1.1"
  if cl=="BE.1.1.1" and "Spike_N460K" in mutations: cl="BQ.1"
  if cl=="BQ.1":
    if "Spike_R346T" in mutations: cl="BQ.1.1"
    elif "Spike_I666V" in mutations: cl="BQ.1.2"
    elif "Spike_E619Q" in mutations: cl="BQ.1.3"
    elif "Spike_R190T" in mutations: cl="BQ.1.4"

  # Offshoots of BA.2.75.*
  if cl=="BA.2.75.2" and "Spike_T604I" in mutations and "Spike_L452R" in mutations: cl="CA.1"
  
  if cl=="BA.2.75.3" and "Spike_F486S" in mutations: cl="BM.1"
  if cl=="BM.1" and "Spike_R346T" in mutations: cl="BM.1.1"
  if cl=="BM.1.1" and "Spike_F490S" in mutations: cl="BM.1.1.1"

  if cl=="BA.2.75.4" and "Spike_R346T" in mutations and "Spike_F486I" in mutations: cl="BR.2"

  if cl=="BA.2.3" and "Spike_E484R" in mutations: cl="BA.2.3.20"
  
  # Recombinant
  if "Spike_G252V" in mutations and ("Spike_F486S" in mutations)+("Spike_F490S" in mutations)+("Spike_R346T" in mutations)>=2: cl="XBB.1"
  
  return cl

# Return classification cl1 if (cl0 is Unassigned) or (cl0 is a prefix of cl1) or (cl1 is a recombinant and cl0 isn't), otherwise return cl0
# So cl0 takes priority if there is a disagreement
def join(cl0,cl1):
  if cl1=="Other": cl1="Unassigned"
  if cl0=="Unassigned": return cl1
  if cl1[0]=="X" and cl0[0]!="X": return cl1
  el0=expandlin(cl0)
  el1=expandlin(cl1)
  if el1[:len(el0)]==el0: return cl1
  return cl0

def classify(mutations,cl="",gisaid=False):
  if gisaid:
    cl=join(cl,classifygisaid.treeclassify(mutations))
    cl=join(cl,manualclassifygisaid(cl,mutations))
  else:
    cl=join(cl,classifycog.treeclassify(mutations))
    cl=join(cl,manualclassifycog(cl,mutations))
  return cl
