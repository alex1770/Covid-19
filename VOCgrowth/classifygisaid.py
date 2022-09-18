lineages = ['B.1.617.2*', 'BA.2.12.1*', 'BA.2.74*', 'BA.2.75*', 'BA.2.76*', 'BA.1*', 'BA.2*', 'BA.3*', 'BA.4*', 'BA.5*', 'Other']

def treeclassify(mutations):
  # Classification run at 2022-07-31.23:09:03
  # Command line: learnclassifier_tree.py -f 2022-06-01 -g -l B.1.617.2*,BA.2.12.1*,BA.2.74*,BA.2.75*,BA.2.76*,BA.1*,BA.2*,BA.3*,BA.4*,BA.5* -m 200 -M 100 -p 30
  # Using input file: metadata_sorted.tsv
  # From: 2022-06-01
  # To: 9999-12-31
  # Mincount: 200
  # Maxleaves: 100
  # Database: GISAID
  # synSNP permissiveness: 1
  # Max N-content: 0.05
  # Number of leaves in decision tree: 30
  if "NSP4_L438F" in mutations:
    if "Spike_S704L" in mutations:
      if "Spike_L452Q" in mutations:
        return "BA.2.12.1*"
      else:
        if "Spike_N440K" in mutations:
          return "BA.2*"
        else:
          return "BA.2.12.1*"
    else:
      if "Spike_Y248N" in mutations:
        return "BA.2.76*"
      else:
        if "Spike_G339H" in mutations:
          return "BA.2.75*"
        else:
          if "Spike_N764K" in mutations:
            if "NS3_P240S" in mutations:
              return "BA.2.74*"
            else:
              if "Spike_F486V" in mutations:
                return "BA.5*"
              else:
                if "Spike_L452Q" in mutations:
                  return "BA.2*"
                else:
                  if "NSP3_G489S" in mutations:
                    if "M_A63T" in mutations:
                      return "BA.2*"
                    else:
                      return "BA.2*"
                  else:
                    return "BA.2*"
          else:
            if "M_A63T" in mutations:
              if "Spike_L452Q" in mutations:
                return "BA.2.12.1*"
              else:
                return "BA.2*"
            else:
              return "BA.2.12.1*"
  else:
    if "NS7b_L11F" in mutations:
      if "N_P151S" in mutations:
        if "Spike_F486V" in mutations:
          return "BA.4*"
        else:
          return "BA.4*"
      else:
        return "BA.5*"
    else:
      if "Spike_F486V" in mutations:
        if "M_D3N" in mutations:
          return "BA.5*"
        else:
          if "N_P151S" in mutations:
            return "BA.4*"
          else:
            if "Spike_H69del" in mutations:
              return "BA.5*"
            else:
              return "BA.5*"
      else:
        if "M_D3N" in mutations:
          if "Spike_L452R" in mutations:
            return "BA.5*"
          else:
            if "Spike_T19I" in mutations:
              if "Spike_R408S" in mutations:
                return "BA.2*"
              else:
                return "BA.5*"
            else:
              return "BA.5*"
        else:
          if "Spike_N856K" in mutations:
            return "BA.1*"
          else:
            if "Spike_S704L" in mutations:
              return "BA.2.12.1*"
            else:
              if "Spike_Q493R" in mutations:
                if "NSP3_T24I" in mutations:
                  return "BA.2*"
                else:
                  return "Other"
              else:
                return "BA.2*"
