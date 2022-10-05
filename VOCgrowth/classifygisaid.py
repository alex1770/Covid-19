lineages = ['BA.5.1.1', 'BA.5.1.2', 'BA.5.3.1', 'BF.11', 'BA.2', 'BA.4.1', 'BA.5.1.5', 'BA.5.2.3', 'BF.14', 'BA.5.9', 'BA.5.2.6', 'BA.5.1.10', 'BE.1', 'BA.5.5', 'BA.5.6', 'BA.5', 'BF.10', 'BA.2.75', 'BF.5', 'BE.1.1', 'BF.7', 'BA.4.6', 'BA.5.2', 'BA.5*', 'Unassigned']
mindate = "2022-08-01"

def treeclassify(mutations):
  # Classification run at 2022-10-05.14:19:25
  # Command line: ./learnclassifier_tree.py -f 2022-08-01 -l BA.5.1.1,BA.5.1.2,BA.5.3.1,BF.11,BA.2,BA.4.1,BA.5.1.5,BA.5.2.3,BF.14,BA.5.9,BA.5.2.6,BA.5.1.10,BE.1,BA.5.5,BA.5.6,BA.5,BF.10,BA.2.75,BF.5,BE.1.1,BF.7,BA.4.6,BA.5.2,BA.5* -p 50 -M 100 -m 50 -g
  # Using input file: metadata_sorted.tsv
  # From: 2022-08-01
  # To: 9999-12-31
  # Mincount: 50
  # Maxleaves: 100
  # Database: GISAID
  # synSNP permissiveness: 1
  # Max N-content: 0.05
  # Number of leaves in decision tree: 50
  if "NSP13_T127N" in mutations:
    if "NSP1_K120N" in mutations:
      return "BA.5.2.3"
    else:
      if "Spike_R346T" in mutations:
        return "BA.5.2.6"
      else:
        if "NSP15_G246S" in mutations:
          return "BA.5*"
        else:
          if "Spike_K444M" in mutations:
            return "BA.5*"
          else:
            return "BA.5.2"
  else:
    if "NS7b_L11F" in mutations:
      if "Spike_R346T" in mutations:
        if "Spike_N658S" in mutations:
          if "N_P151S" in mutations:
            return "BA.4.6"
          else:
            return "BA.5"
        else:
          return "BA.4.1"
      else:
        if "Spike_V3G" in mutations:
          return "BA.4.1"
        else:
          return "Unassigned"
    else:
      if "NSP13_M233I" in mutations:
        if "NSP6_L260F" in mutations:
          if "Spike_F486V" in mutations:
            return "BE.1.1"
          else:
            return "BA.2"
        else:
          if "Spike_R346T" in mutations:
            return "BA.5*"
          else:
            if "N_E136D" in mutations:
              return "BE.1"
            else:
              return "BA.5*"
      else:
        if "Spike_A1020S" in mutations:
          if "NS7a_H47Y" in mutations:
            return "BF.5"
          else:
            return "BA.5*"
        else:
          if "Spike_T76I" in mutations:
            if "Spike_N450D" in mutations:
              return "BA.5*"
            else:
              return "BA.5.5"
          else:
            if "NSP13_A389V" in mutations:
              if "Spike_K444T" in mutations:
                return "BA.5*"
              else:
                if "NSP8_M174V" in mutations:
                  return "BA.5*"
                else:
                  return "BA.5.6"
            else:
              if "NSP4_L438F" in mutations:
                if "NSP8_N118S" in mutations:
                  return "BA.2.75"
                else:
                  return "Unassigned"
              else:
                if "Spike_R346T" in mutations:
                  if "N_S33del" in mutations:
                    if "N_P151S" in mutations:
                      return "Unassigned"
                    else:
                      return "BA.5*"
                  else:
                    return "BF.7"
                else:
                  if "N_P13F" in mutations:
                    if "NSP13_A368V" in mutations:
                      return "BF.10"
                    else:
                      return "BA.5*"
                  else:
                    if "NS3_L94I" in mutations:
                      return "BA.5.1.1"
                    else:
                      if "NSP12_Y521C" in mutations:
                        if "Spike_R346I" in mutations:
                          return "BA.5.9"
                        else:
                          return "BA.5.2"
                      else:
                        if "NSP6_K270R" in mutations:
                          if "NSP2_G339S" in mutations:
                            return "BA.5*"
                          else:
                            return "BA.5.1.10"
                        else:
                          if "NS6_V5I" in mutations:
                            return "BA.5.1.2"
                          else:
                            if "N_E136D" in mutations:
                              if "NSP13_Y324C" in mutations:
                                return "BA.5*"
                              else:
                                if "NSP3_T64I" in mutations:
                                  return "BA.5*"
                                else:
                                  return "BA.5.3.1"
                            else:
                              if "Spike_P1263Q" in mutations:
                                return "BA.5.1.5"
                              else:
                                if "Spike_N450D" in mutations:
                                  return "BF.14"
                                else:
                                  if "N_G30del" in mutations:
                                    return "BA.5*"
                                  else:
                                    if "M_D3N" in mutations:
                                      if "NSP3_A1819T" in mutations:
                                        return "BA.5"
                                      else:
                                        if "NSP3_P200S" in mutations:
                                          return "BA.5"
                                        else:
                                          if "Spike_R346I" in mutations:
                                            return "BA.5.9"
                                          else:
                                            if "Spike_D138N" in mutations:
                                              return "BA.5.2"
                                            else:
                                              if "NSP14_A119V" in mutations:
                                                return "BA.5"
                                              else:
                                                return "BA.5*"
                                    else:
                                      if "Spike_F486V" in mutations:
                                        if "N_P151S" in mutations:
                                          return "BA.4.1"
                                        else:
                                          return "BA.5*"
                                      else:
                                        return "Unassigned"
