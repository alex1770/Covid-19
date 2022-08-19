lineages = ['BA.5.1', 'BA.5.2', 'BA.5.2.1', 'BA.4', 'BA.4.1', 'BE.1', 'BE.1.1', 'BA.4.6', 'BA.5', 'BA.5.2.3', 'BA.2.12.1', 'BF.5', 'BF.1', 'BA.5.3.3', 'BA.5.5', 'BA.5.8', 'BF.6', 'BA.2', 'BF.4', 'BA.5.6', 'BA.5.3.1', 'BA.5.1.5', 'BA.5.1.2', 'BF.11', 'BA.4.7', 'BA.5.1.3', 'BF.10', 'BF.7', 'BA.5.1.4', 'BA.5.1.1', 'BA.2.18', 'BA.4.1.1', 'BA.2.76', 'BA.5.2.2', 'BE.3', 'BA.5.3', 'BA.4.1.8', 'BA.2.75', 'BF.2', 'BA.2.56', 'BA.2.9', 'BA.5.9', 'BA.2.23', 'BA.2.38', 'BF.3', 'BA.2.36', 'BA.5.3.2', 'BA.4.5', 'BA.2.3', 'BA.4.1.6', 'Other']
mindate = "2022-07-01"

def treeclassify(mutations):
  # Classification run at 2022-08-19.00:25:10
  # Command line: learnclassifier_tree.py -f 2022-07-01 -n50 -p 50 -M 100 -m 10
  # Using input file: cog_metadata_sorted.csv
  # From: 2022-07-01
  # To: 9999-12-31
  # Mincount: 10
  # Maxleaves: 100
  # Database: COG-UK
  # synSNP permissiveness: 1
  # Max N-content: 0.05
  # Number of leaves in decision tree: 50
  if "synSNP:A28330G" in mutations:
    if "orf1ab:T5451N" in mutations:
      if "orf1ab:K120N" in mutations:
        return "BA.5.2.3"
      else:
        if "orf1ab:G6698S" in mutations:
          return "BA.5.2.2"
        else:
          return "BA.5.2"
    else:
      if "S:A1020S" in mutations:
        if "ORF7a:H47Y" in mutations:
          return "BF.5"
        else:
          return "BF.3"
      else:
        if "orf1ab:V7086F" in mutations:
          return "BF.1"
        else:
          if "S:G181A" in mutations:
            return "BF.6"
          else:
            if "S:R346T" in mutations:
              if "N:S33F" in mutations:
                return "BF.7"
              else:
                return "BF.11"
            else:
              if "S:T259A" in mutations:
                return "BF.4"
              else:
                if "orf1ab:A5692V" in mutations:
                  return "BF.10"
                else:
                  if "orf1ab:Y4913C" in mutations:
                    return "BA.5.2"
                  else:
                    if "orf1ab:V4712L" in mutations:
                      return "BF.2"
                    else:
                      if "N:S33F" in mutations:
                        return "BA.5.2.1"
                      else:
                        return "BA.5.2.1"
  else:
    if "ORF10:L37F" in mutations:
      if "ORF6:V5I" in mutations:
        return "BA.5.1.2"
      else:
        if "synSNP:T28510C" in mutations:
          return "BA.5.1.5"
        else:
          if "S:V289I" in mutations:
            return "BA.5.1.3"
          else:
            if "orf1ab:V1117I" in mutations:
              return "BA.5.1.4"
            else:
              if "ORF3a:L94I" in mutations:
                return "BA.5.1.1"
              else:
                if "ORF6:D61L" in mutations:
                  return "BA.5"
                else:
                  return "BA.5.1"
    else:
      if "N:P151S" in mutations:
        if "S:V3G" in mutations:
          if "ORF3a:I37T" in mutations:
            return "BA.4.1.1"
          else:
            if "S:R346T" in mutations:
              return "BA.4.1.8"
            else:
              if "orf1ab:T4355I" in mutations:
                return "BA.4.1.6"
              else:
                return "BA.4.1"
        else:
          if "S:R346T" in mutations:
            return "BA.4.6"
          else:
            if "ORF3a:I47V" in mutations:
              return "BA.4.7"
            else:
              if "M:D3N" in mutations:
                return "BE.3"
              else:
                if "S:Q1201L" in mutations:
                  return "BA.4.5"
                else:
                  return "BA.4"
      else:
        if "orf1ab:M5557I" in mutations:
          if "orf1ab:L3829F" in mutations:
            return "BE.1.1"
          else:
            return "BE.1"
        else:
          if "S:F486V" in mutations:
            if "orf1ab:Q556K" in mutations:
              if "orf1ab:S2193F" in mutations:
                return "BA.5.3.3"
              else:
                if "N:E136D" in mutations:
                  return "BA.5.3.1"
                else:
                  return "BA.5.3"
            else:
              if "S:T76I" in mutations:
                return "BA.5.5"
              else:
                if "S:P1162L" in mutations:
                  return "BA.5.8"
                else:
                  if "orf1ab:A5713V" in mutations:
                    return "BA.5.6"
                  else:
                    if "S:R346I" in mutations:
                      return "BA.5.9"
                    else:
                      if "ORF6:D61L" in mutations:
                        return "BA.2"
                      else:
                        if "S:P1263Q" in mutations:
                          return "BA.5.1.5"
                        else:
                          return "BA.5"
          else:
            if "S:S704L" in mutations:
              return "BA.2.12.1"
            else:
              if "S:K417T" in mutations:
                return "BA.2.18"
              else:
                if "S:Y248N" in mutations:
                  return "BA.2.76"
                else:
                  if "S:W152R" in mutations:
                    return "BA.2.75"
                  else:
                    if "ORF3a:H78Y" in mutations:
                      return "BA.2.9"
                    else:
                      if "S:L452M" in mutations:
                        return "BA.2.56"
                      else:
                        return "BA.2"
