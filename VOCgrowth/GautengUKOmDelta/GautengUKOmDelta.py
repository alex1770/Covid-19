# Sources (as at 2021-12-06):
# https://sacoronavirus.co.za/latest-vaccine-statistics/
# https://coronavirus.data.gov.uk/details/vaccinations?areaType=overview&areaName=United%20Kingdom
# https://www.gov.uk/government/publications/covid-19-vaccine-weekly-surveillance-reports
# http://sonorouschocolate.com/covid19/index.php?title=Early_Omicron_Growth_Estimate
#
# See https://docs.google.com/spreadsheets/d/1wvQdozRGTpCvU44IQiT8D5z6dQVBR-9xTHlZ5vhFzW0/edit#gid=812895892
# for a spreadsheet illustration of this program

from math import exp,log


# Input is vaccine rate array, attack rate. VR[i] = proportion of population who have had exactly i doses
# Assume vaccination status is independent of whether infected (yes, I know)
def IndImStatus(VR,AR):
  return [VR[0]*(1-AR),VR[0]*AR,
          VR[1]*(1-AR),VR[1]*AR,
          VR[2]*(1-AR),VR[2]*AR,
          VR[3]]

def getdoublingtime(
    # Vaccination rates
    # 0 dose, 1 dose, 2 dose, 3 dose
    VR_Gauteng=[0.430, 0.364, 0.206, 0.],
    VR_UK=[0.238, 0.068, 0.392, 0.302],
    
    # Attack rates
    AR_Gauteng=0.95,
    AR_UK=0.45,
    
    # Generation times
    T_Del=5,
    T_Om=5,

    # Growth rates
    g_Del_Gauteng=-0.043,
    g_Om_Gauteng=0.232,
    g_Del_UK=0.005,

    # Immunity coefficients vs new symptomatic infection
    # Map from immunity status as defined by
    # [Naive, pre-infected only, 1 dose, pre-infected + 1 dose, 2 dose, pre-infected + 2 dose, 3 dose]
    # Can also specify Im_Om as a logistic shift from Im_Del
    Im_Del=[0, 0.75, 0.30, 0.80, 0.60, 0.92, 0.92],
    Im_Om=0.1#[0, 0.15, 0.05, 0.40, 0.30, 0.70, 0.80]
  ):
  
  assert abs(sum(VR_Gauteng)-1)<1e-6 and abs(sum(VR_UK)-1)<1e-6
  ImStatus_Gauteng=IndImStatus(VR_Gauteng,AR_Gauteng)
  ImStatus_UK=IndImStatus(VR_UK,AR_UK)

  if type(Im_Om)!=list:
    Im_Om=[Im_Om*p/(Im_Om*p+1-p) for p in Im_Del]
  
  RR_Gauteng_Del=1-sum(x*y for (x,y) in zip(ImStatus_Gauteng,Im_Del))
  RR_Gauteng_Om=1-sum(x*y for (x,y) in zip(ImStatus_Gauteng,Im_Om))
  RR_UK_Del=1-sum(x*y for (x,y) in zip(ImStatus_UK,Im_Del))
  RR_UK_Om=1-sum(x*y for (x,y) in zip(ImStatus_UK,Im_Om))
  
  RRR_Gauteng=RR_Gauteng_Om/RR_Gauteng_Del
  Rrel_obs_Gauteng=exp(T_Om*g_Om_Gauteng-T_Del*g_Del_Gauteng)
  RR0=Rrel_obs_Gauteng/RRR_Gauteng
  RRR_UK=RR_UK_Om/RR_UK_Del
  Rrel_pred_UK=RR0*RRR_UK
  # Rrel_pred_UK=exp(T_Om*g_Om_UK-T_Del*g_Del_UK)
  g_Om_UK=(log(Rrel_pred_UK)+T_Del*g_Del_UK)/T_Om
  return log(2)/g_Om_UK,RR0

