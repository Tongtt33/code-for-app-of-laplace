library(VGAM)
library(stats4)
library(splines)
library(gmp)
library(hms)

shinyServer(function(input, output) {
  
  output$values <- renderTable({
    
    
    
    
    findnmle<-function(f,c){
      
      vec<-matrix(1:181,ncol=2,nrow= 181)
      
      for(m in 20:200){
        #c <- 0.99
       # f <- 0.1
        
        L<-(-f)*sqrt(2)
        
        U<-f*sqrt(2)
        
        n <- 2*m-1
        pdf_l <-  function(x){
          con <- factorialZ(n)/(factorialZ(m-1)*factorialZ(n-m))
          
          val <- (con)*dlaplace(x,0,1)*(plaplace(x,0,1)^(m-1))*((1-plaplace(x,0,1))^(m-1))
          value <- as.double(val)
        }
        
        cdf_l <- function(q){
          
          value <- integrate(pdf_l, lower = -Inf, upper = q)$value
          
          return(value)
          
        }
        
        vec[m-19,1]<- n
        
        vec[m-19,2]<- cdf_l(U)-(1+c)/2
        
      }
      
      vec1<- vec[order(vec[,2]),]
      
      vec2 <- vec1[which( vec1[,2]>=0 ),]
      
      return(vec2[1,1])
      
    }
    
    findnme<-function(f,c){
      
      vec_1<-matrix(1:40,ncol=2,nrow=40  )
      vec_2 <- matrix(1:11, ncol=2,nrow=11)
       #c <- 0.95
       #f <- 0.1
      #n <- 30
      for(k in 1:40){
        n <- k*10
        L<-(-f)*sqrt(2)*n
        
        U<-f*sqrt(2)*n
        pdf_W <- function(w){
          val <- dgamma(w^2,shape = n,rate = 1/2)*2*w
        }
        
        pdf_U <-  function(u){
          valfunc <- function(w,u){
            
            inte <- dnorm(u/w,0,1)*pdf_W(w)/w
          }
          value <- integrate(Vectorize(valfunc),lower = 0,upper = 150,u=u)$value
          final_value <- value
          return(final_value)
        }
        
        cdf_U <- function(q){
          val <- integrate(Vectorize(pdf_U),lower = -Inf,upper = q)$value
        }
        
        vec_1[k,1]<- k
        
        vec_1[k,2]<- cdf_U(U)-cdf_U(L)-c
        
      }
      
      vec_11<- vec_1[order(vec_1[,2]),]
      
      vec_12 <- vec_11[which( vec_11[,2]>=0 ),]
     
######
      for(j in 0:10){
        n <- (vec_12[1]-1)*10+j
        L<-(-f)*sqrt(2)*n
        
        U<-f*sqrt(2)*n
        pdf_W <- function(w){
          val <- dgamma(w^2,shape = n,rate = 1/2)*2*w
        }
        
        pdf_U <-  function(u){
          valfunc <- function(w,u){
            
            inte <- dnorm(u/w,0,1)*pdf_W(w)/w
          }
          value <- integrate(Vectorize(valfunc),lower = 0,upper = 150,u=u)$value
          final_value <- value
          return(final_value)
        }
        
        cdf_U <- function(q){
          val <- integrate(Vectorize(pdf_U),lower = -Inf,upper = q)$value
        }
        
        vec_2[j+1,1]<- n
        
        vec_2[j+1,2]<- cdf_U(U)-cdf_U(L)-c
        
      }
      
      vec_21<- vec_2[order(vec_2[,2]),]
      
      vec_22 <- vec_21[which( vec_21[,2]>=0 ),]
      
      
      
      return(vec_22[1,1])
      
    }
    
    timestart<-Sys.time()
    n_mle<-findnmle(input$f,input$c)
    n_me<-findnme(input$f,input$c)
    timeend<-Sys.time()
    duration<-floor(difftime(timeend,timestart,units="secs"))
    runningtime<-as.character(as_hms(duration))
    Outcome<- structure(c("Precision (f)","Confidence Level (c)","Sample Size when using sample median (n_mle)","Sample Size when using sample mean (n_me)","Running time",
                          input$f, input$c, n_mle, n_me, runningtime ),.Dim=c(5L,2L))
    colnames(Outcome)<-c("Name","Value")
    return (Outcome)
    
  })
})

