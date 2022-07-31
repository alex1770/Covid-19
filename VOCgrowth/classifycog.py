lineages = ['B.1.617.2*', 'BA.2.12.1*', 'BA.2.74*', 'BA.2.75*', 'BA.2.76*', 'BA.1*', 'BA.2*', 'BA.3*', 'BA.4*', 'BA.5*', 'Other']

def treeclassify(mutations):
  # Classification run at 2022-07-31.23:01:18
  # Command line: learnclassifier_tree.py -f 2022-06-01 -l B.1.617.2*,BA.2.12.1*,BA.2.74*,BA.2.75*,BA.2.76*,BA.1*,BA.2*,BA.3*,BA.4*,BA.5* -m 50 -M 100 -p 30
  # Using input file: cog_metadata_sorted.csv
  # From: 2022-06-01
  # To: 9999-12-31
  # Mincount: 50
  # Maxleaves: 100
  # Database: COG-UK
  # synSNP permissiveness: 1
  # Max N-content: 0.05
  # Number of leaves in decision tree: 30
  if "M:D3N" in mutations:
    if "orf1ab:L3201F" in mutations:
      return "BA.2*"
    else:
      if "ORF6:D61L" in mutations:
        return "BA.5*"
      else:
        return "BA.5*"
  else:
    if "N:P151S" in mutations:
      if "ORF6:D61L" in mutations:
        if "orf1ab:T6564I" in mutations:
          if "orf1ab:L3201F" in mutations:
            return "BA.2*"
          else:
            if "S:S371F" in mutations:
              if "orf1ab:G1307S" in mutations:
                return "BA.4*"
              else:
                return "BA.4*"
            else:
              return "BA.4*"
        else:
          if "S:R346T" in mutations:
            if "S:N658S" in mutations:
              return "BA.4*"
            else:
              return "BA.2*"
          else:
            if "orf1ab:L3201F" in mutations:
              return "BA.2*"
            else:
              return "BA.4*"
      else:
        if "E:V62F" in mutations:
          return "BA.4*"
        else:
          return "BA.4*"
    else:
      if "S:L452Q" in mutations:
        if "S:S704L" in mutations:
          return "BA.2.12.1*"
        else:
          return "BA.2*"
      else:
        if "orf1ab:L3201F" in mutations:
          if "S:Y248N" in mutations:
            return "BA.2.76*"
          else:
            if "S:G446S" in mutations:
              return "BA.2.75*"
            else:
              if "S:S704L" in mutations:
                if "S:R408S" in mutations:
                  return "BA.2*"
                else:
                  return "BA.2.12.1*"
              else:
                if "S:R346T" in mutations:
                  if "S:L452M" in mutations:
                    return "BA.2.74*"
                  else:
                    return "BA.2*"
                else:
                  if "S:N764K" in mutations:
                    return "BA.2*"
                  else:
                    return "BA.2*"
        else:
          if "S:F486V" in mutations:
            if "ORF6:D61L" in mutations:
              if "S:S371F" in mutations:
                if "orf1ab:T6564I" in mutations:
                  return "BA.5*"
                else:
                  return "BA.2*"
              else:
                return "BA.2*"
            else:
              return "BA.5*"
          else:
            if "orf1ab:G1307S" in mutations:
              return "BA.2*"
            else:
              if "S:A67V" in mutations:
                return "BA.1*"
              else:
                return "Other"
