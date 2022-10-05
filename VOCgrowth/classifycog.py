lineages = ['BA.5.1.2', 'BA.5.1.12', 'BL.1', 'BA.4.7', 'BA.5.1.10', 'BA.5.6', 'BE.1', 'BA.2.75.3', 'BF.10', 'BA.5.3.1', 'BA.5.1.5', 'BA.5.2.3', 'BA.2.75.5', 'BA.5', 'BA.2.75.2', 'BA.5.9', 'BF.5', 'BF.11', 'BA.5.2.6', 'BF.7', 'BA.4.6', 'BE.1.1', 'BA.5.2.1', 'BA.5.1', 'BA.5.2', 'Other']
mindate = "2022-08-01"

def treeclassify(mutations):
  # Classification run at 2022-10-05.10:17:09
  # Command line: ./learnclassifier_tree.py -f 2022-08-01 -l BA.5.1.2,BA.5.1.12,BL.1,BA.4.7,BA.5.1.10,BA.5.6,BE.1,BA.2.75.3,BF.10,BA.5.3.1,BA.5.1.5,BA.5.2.3,BA.2.75.5,BA.5,BA.2.75.2,BA.5.9,BF.5,BF.11,BA.5.2.6,BF.7,BA.4.6,BE.1.1,BA.5.2.1,BA.5.1,BA.5.2 -p 50 -M 100 -m 10
  # Using input file: cog_metadata_sorted.csv
  # From: 2022-08-01
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
        if "S:R346T" in mutations:
          return "BA.5.2.6"
        else:
          if "S:K444M" in mutations:
            return "Other"
          else:
            if "orf1ab:G6698S" in mutations:
              return "Other"
            else:
              return "BA.5.2"
    else:
      if "S:R346T" in mutations:
        if "N:S33F" in mutations:
          return "BF.7"
        else:
          return "BF.11"
      else:
        if "ORF7a:H47Y" in mutations:
          return "BF.5"
        else:
          if "N:P13F" in mutations:
            if "orf1ab:A5692V" in mutations:
              return "BF.10"
            else:
              return "BA.5.2.1"
          else:
            if "orf1ab:Y4913C" in mutations:
              return "BA.5.2"
            else:
              if "S:G181A" in mutations:
                return "Other"
              else:
                if "orf1ab:V7086F" in mutations:
                  return "Other"
                else:
                  if "S:T259A" in mutations:
                    return "Other"
                  else:
                    if "S:T547I" in mutations:
                      return "Other"
                    else:
                      if "S:Y248H" in mutations:
                        return "Other"
                      else:
                        if "N:S33F" in mutations:
                          return "BA.5.2.1"
                        else:
                          if "orf1ab:V4712L" in mutations:
                            return "Other"
                          else:
                            if "orf1ab:V3858L" in mutations:
                              return "Other"
                            else:
                              if "S:N450D" in mutations:
                                return "Other"
                              else:
                                return "BA.5.2.1"
  else:
    if "ORF10:L37F" in mutations:
      if "synSNP:T28510C" in mutations:
        return "BA.5.1.5"
      else:
        if "ORF6:V5I" in mutations:
          return "BA.5.1.2"
        else:
          if "orf1ab:K3839R" in mutations:
            return "BA.5.1.10"
          else:
            if "S:V445A" in mutations:
              return "BA.5.1.12"
            else:
              if "S:V289I" in mutations:
                return "Other"
              else:
                if "ORF3a:L94I" in mutations:
                  return "Other"
                else:
                  if "orf1ab:V1117I" in mutations:
                    return "Other"
                  else:
                    if "orf1ab:M4855V" in mutations:
                      return "Other"
                    else:
                      return "BA.5.1"
    else:
      if "orf1ab:M5557I" in mutations:
        if "orf1ab:L3829F" in mutations:
          return "BE.1.1"
        else:
          if "S:R346T" in mutations:
            return "Other"
          else:
            return "BE.1"
      else:
        if "S:N658S" in mutations:
          if "S:R346T" in mutations:
            return "BA.4.6"
          else:
            if "ORF3a:I47V" in mutations:
              return "BA.4.7"
            else:
              return "Other"
        else:
          if "M:D3N" in mutations:
            if "orf1ab:A5713V" in mutations:
              return "BA.5.6"
            else:
              if "S:R346I" in mutations:
                return "BA.5.9"
              else:
                if "N:E136D" in mutations:
                  if "orf1ab:Y5648C" in mutations:
                    return "Other"
                  else:
                    return "BA.5.3.1"
                else:
                  if "orf1ab:Q556K" in mutations:
                    return "Other"
                  else:
                    if "S:T76I" in mutations:
                      return "Other"
                    else:
                      if "S:P1162L" in mutations:
                        return "Other"
                      else:
                        if "orf1ab:Q3922R" in mutations:
                          return "Other"
                        else:
                          return "BA.5"
          else:
            if "S:K147E" in mutations:
              if "S:F486S" in mutations:
                if "orf1ab:V6107I" in mutations:
                  return "BA.2.75.3"
                else:
                  return "BA.2.75.2"
              else:
                if "S:K356T" in mutations:
                  return "BA.2.75.5"
                else:
                  return "Other"
            else:
              return "Other"
