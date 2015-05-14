##Functions
  #function used to extract specific range of data
  checkThreshold <- function(i, cond1, cond2){
    if(!is.na(QP$Rawdat.delta_P[i]) && QP$Rawdat.delta_P[i]>=cond1 && QP$Rawdat.delta_P[i]<cond2){
      newVar <- QP$Rawdat.delta_P[i]
      return(newVar)
    }
  }
  
  #function used to extract specific range of data matching the row of check Threshold()
  checkThreshold2 <- function(i, cond1, cond2){
    if(!is.na(QP$Rawdat.delta_P[i]) && QP$Rawdat.delta_P[i]>=cond1 && QP$Rawdat.delta_P[i]<cond2){
      newVar <- QP$Rawdat.ome_Q[i]
      return(newVar)
    }
  }
  
  #function used to obtain calibration coefficient C using bisection method for subset of data
  obtainCal <- function(dataSet, iter){
    if(length(dataSet)>0){
      c_low <- -10
      c_high <- 10
      for(i in 1:iter){
        ##Calculates Flow rate using high and low calibration coefficient
        fixedQ <- mean(dataSet[,c(2)])
        qlow <- constA*c_low*dataSet[,c(1)]
        qhigh <- constA*c_high*dataSet[,c(1)]
        qavg <- mean(rbind(qlow,qhigh))
        cavg <- (c_low + c_high)/2
        #print(cavg)
        #print(qavg)
        #print(fixedQ)
        
        if((fixedQ-qavg)<0){
          c_high=cavg
        }
        if((fixedQ-qavg)>0){
          c_low=cavg
        }
      }
    }
    else{cavg<-0}
    cavg
  }
  
  #function calculating percent difference between calibrated value and original value
  percent_difference <- function(ori_dat1, ori_dat2, cal_file, del_1, del_2){
    range_0 = paste(as.character(del_1), "-", as.character(del_2))
    if(!is.null(ori_dat1) && !is.null(ori_dat2)){  
      per_diff<-c()
      per_average<-c()
      per_max<-c()
      per_min<-c()
      low_limit <- c()
      high_limit <- c()
      range <- c()
      for (i in 1:length(cal_file)){
        cal_value <- cal_file[i]
        cal_flow <- ori_dat1*cal_value*constA
        low_limit_0 <- del_1*cal_value*constA
        high_limit_0 <- del_2*cal_value*constA
        low_limit <- c(low_limit, low_limit_0)
        high_limit <- c(high_limit, high_limit_0)
        
        for (k in 1:length(cal_flow)){
          per_diff <- c(per_diff, ((cal_flow[k]-ori_dat2[k])/ori_dat2[k])*100)
        }
        per_average <- c(per_average, mean(per_diff))
        per_max <- c(per_max, max(per_diff))
        per_min <- c(per_min, min(per_diff))
        range <- c(range, range_0)
      }
      result <- data.frame(range, cal_file, per_average, per_max, per_min, low_limit, high_limit)
    }
    else
      result <- c("NA","NA", "NA", "NA", "NA","NA","NA")
  }
  
##Global Variables
  #variable Required
  cal_1 <- c()
  cal_2 <- c()
  cal_3 <- c()
  cal_4 <- c()
  cal_5 <- c()
  cal_6 <- c()
  cal_7 <- c()
  cal_8 <- c()
  cal_9 <- c()
  cal_10 <- c()
  cal_11 <- c()
  cal_12 <- c()
  cal_13 <- c()
  cal_14 <- c()
  cal_15 <- c()
  cal_16 <- c()
  cal_17 <- c()
  cal_18 <- c()
  cal_19 <- c()
  cal_20 <- c()
  
  Dia <- 0.0147 #pipe diameter in m
  dia <- 0.0104 #Orifice diameter in m 
  beta <- dia/Dia
  temp <- 20+274.15 #Gas Temperature Degree Kelvin
  Pressure <- 101.325 #KPa Atmospheric
  in_P <- 7*0.248126 #KPa Gas Pressure
  Gas_R <- 0.5182 #Universal Gas Constant
  row <- (in_P+Pressure)/(temp*Gas_R)
  A_1 <- 1/(sqrt(1-beta^4))
  A_2 <- (pi/4)*(dia^2)
  A_3 <- sqrt(2*row)
  A_4 <- (60*35.31467/row)
  constA <- A_1*A_2*A_3*A_4 #geometric values - orifice
  
  Rawdatfile <- list.files("cal_air", full.names=TRUE) #"specdata" working directory within the folder


## Main Program

for(i in 1:length(Rawdatfile)){ #use length(Rawdatfile) later on
  #initialize all local variable
  Q1 <- c()
  Q2 <- c()
  Q3 <- c()
  Q4 <- c()
  Q5 <- c()
  Q6 <- c()
  Q7 <- c()
  Q8 <- c()
  Q9 <- c()
  Q10 <- c()
  Q11 <- c()
  Q12 <- c()
  Q13 <- c()
  Q14 <- c()
  Q15 <- c()
  Q16 <- c()
  Q17 <- c()
  Q18 <- c()
  Q19 <- c()
  Q20 <- c()
  Q21 <- c()
  Q22 <- c()
  Q23 <- c()
  Q24 <- c()
  Q25 <- c()
  Q26 <- c()
  Q27 <- c()
  Q28 <- c()
  Q29 <- c()
  Q30 <- c()
  Q31 <- c()
  Q32 <- c()
  Q33 <- c()
  Q34 <- c()
  Q35 <- c()
  Q36 <- c()
  Q37 <- c()
  Q38 <- c()
  Q39 <- c()
  Q40 <- c()
  
  Rawdat = read.csv(Rawdatfile[i]) ##Read individual files from folder
  QP = data.frame(Rawdat$delta_P, Rawdat$ome_Q)  #Extract Specific column from matrix
  
  #Loop through values to categorize range in different subset
  for(i in 1:length(QP$Rawdat.delta_P)){
    Q1 <- c(Q1, checkThreshold(i, 0, 20))
    Q2 <- c(Q2, checkThreshold2(i, 0, 20))
    
    Q3 <- c(Q3, checkThreshold(i, 20,40))
    Q4 <- c(Q4, checkThreshold2(i,20,40))
    
    Q5 <- c(Q5, checkThreshold(i, 40,60))
    Q6 <- c(Q6, checkThreshold2(i,40,60)) 
    
    Q7 <- c(Q7, checkThreshold(i, 60,80))
    Q8 <- c(Q8, checkThreshold2(i,60,80))
    
    Q9 <- c(Q9, checkThreshold(i, 80,100))
    Q10 <- c(Q10, checkThreshold2(i,80,100))
    
    Q11 <- c(Q11, checkThreshold(i, 100,120))
    Q12 <- c(Q12, checkThreshold2(i,100,120))
    
    Q13 <- c(Q13, checkThreshold(i, 120,140))
    Q14 <- c(Q14, checkThreshold2(i,120,140))
    
    Q15 <- c(Q15, checkThreshold(i, 140,160))
    Q16 <- c(Q16, checkThreshold2(i,140,160))
    
    Q17 <- c(Q17, checkThreshold(i, 160,180))
    Q18 <- c(Q18, checkThreshold2(i,160,180))
    
    Q19 <- c(Q19, checkThreshold(i, 180,200))
    Q20 <- c(Q20, checkThreshold2(i,180,200))
    
    Q21 <- c(Q21, checkThreshold(i, 200,220))
    Q22 <- c(Q22, checkThreshold2(i,200,220))
    
    Q23 <- c(Q23, checkThreshold(i, 220,240))
    Q24 <- c(Q24, checkThreshold2(i,220,240))
    
    Q25 <- c(Q25, checkThreshold(i, 240,260))
    Q26 <- c(Q26, checkThreshold2(i,240,260))
    
    Q27 <- c(Q27, checkThreshold(i, 260,280))
    Q28 <- c(Q28, checkThreshold2(i,260,280))
    
    Q29 <- c(Q29, checkThreshold(i, 280,300))
    Q30 <- c(Q30, checkThreshold2(i,280,300))
    
    Q31 <- c(Q31, checkThreshold(i, 300,320))
    Q32 <- c(Q32, checkThreshold2(i,300,320))
    
    Q33 <- c(Q33, checkThreshold(i, 320,340))
    Q34 <- c(Q34, checkThreshold2(i,320,340))
    
    Q35 <- c(Q35, checkThreshold(i, 340,360))
    Q36 <- c(Q36, checkThreshold2(i,340,360))
    
    Q37 <- c(Q37, checkThreshold(i, 360,380))
    Q38 <- c(Q38, checkThreshold2(i,360,380))
    
    Q39 <- c(Q39, checkThreshold(i, 380,400))
    Q40 <- c(Q40, checkThreshold2(i,380,400))
  }
  
  dat00_02 <- data.frame(Q1, Q2)
  dat02_04 <- data.frame(Q3, Q4)
  dat04_06 <- data.frame(Q5, Q6)
  dat06_08 <- data.frame(Q7, Q8)
  dat08_10 <- data.frame(Q9, Q10)
  dat10_12 <- data.frame(Q11, Q12)
  dat12_14 <- data.frame(Q13, Q14)
  dat14_16 <- data.frame(Q15, Q16)
  dat16_18 <- data.frame(Q17, Q18)
  dat18_20 <- data.frame(Q19, Q20)
  dat20_22 <- data.frame(Q21, Q22)
  dat22_24 <- data.frame(Q23, Q24)
  dat24_26 <- data.frame(Q25, Q26)
  dat26_28 <- data.frame(Q27, Q28)
  dat28_30 <- data.frame(Q29, Q30)
  dat30_32 <- data.frame(Q31, Q32)
  dat32_34 <- data.frame(Q33, Q34)
  dat34_36 <- data.frame(Q35, Q36)
  dat36_38 <- data.frame(Q37, Q38)
  dat38_40 <- data.frame(Q39, Q40)
  
  cal_1 <- c(cal_1, obtainCal(dat00_02, 100))
  cal_2 <- c(cal_2, obtainCal(dat02_04, 100))
  cal_3 <- c(cal_3, obtainCal(dat04_06, 100))
  cal_4 <- c(cal_4, obtainCal(dat06_08, 100))
  cal_5 <- c(cal_5, obtainCal(dat08_10, 100))
  cal_6 <- c(cal_6, obtainCal(dat10_12, 100))
  cal_7 <- c(cal_7, obtainCal(dat12_14, 100))
  cal_8 <- c(cal_8, obtainCal(dat14_16, 100))
  cal_9 <- c(cal_9, obtainCal(dat16_18, 100))
  cal_10 <- c(cal_10, obtainCal(dat18_20, 100))
  cal_11 <- c(cal_11, obtainCal(dat20_22, 100))
  cal_12 <- c(cal_12, obtainCal(dat22_24, 100))
  cal_13 <- c(cal_13, obtainCal(dat24_26, 100))
  cal_14 <- c(cal_14, obtainCal(dat26_28, 100))
  cal_15 <- c(cal_15, obtainCal(dat28_30, 100))
  cal_16 <- c(cal_16, obtainCal(dat30_32, 100))
  cal_17 <- c(cal_17, obtainCal(dat32_34, 100))
  cal_18 <- c(cal_18, obtainCal(dat34_36, 100))
  cal_19 <- c(cal_19, obtainCal(dat36_38, 100))
  cal_20 <- c(cal_20, obtainCal(dat38_40, 100))
}

  per1 <- percent_difference(Q1, Q2, cal_1, 0, 20)
  per2 <- percent_difference(Q3, Q4, cal_2, 20, 40)
  per3 <- percent_difference(Q5, Q6, cal_3, 40, 60)
  per4 <- percent_difference(Q7, Q8, cal_4, 60, 80)
  per5 <- percent_difference(Q9, Q10, cal_5, 80, 100)
  per6 <- percent_difference(Q11, Q12, cal_6, 100, 120)
  per7 <- percent_difference(Q13, Q14, cal_7, 120, 140)
  per8 <- percent_difference(Q15, Q16, cal_8, 140, 160)
  per9 <-  percent_difference(Q17, Q18, cal_9, 160, 180)
  per10 <- percent_difference(Q19, Q20, cal_10, 180, 200)
  per11 <- percent_difference(Q21, Q22, cal_11, 200, 220)
  per12 <- percent_difference(Q23, Q24, cal_12, 220, 240)
  per13 <-  percent_difference(Q25, Q26, cal_13, 240, 260)
  per14 <-  percent_difference(Q27, Q28, cal_14, 260, 280)
  per15 <-  percent_difference(Q29, Q30, cal_15, 280, 300)
  per16 <-  percent_difference(Q31, Q32, cal_16, 300, 320)
  per17 <-  percent_difference(Q33, Q34, cal_17, 320, 340)
  per18 <-  percent_difference(Q35, Q36, cal_18, 340, 360)
  per19 <- percent_difference(Q37, Q38, cal_19, 360, 380)
  per20 <-  percent_difference(Q39, Q40, cal_20, 380, 400)
  
  final_result <- rbind(per1,per2,per3,per4,per5,per6,per7,per8,per9,per10,per11,per12,per13,per14,per15,per16,per17,per18,per19,per20)
  print(final_result)
  
  ##write.csv(MyData, file = "MyData.csv")
  