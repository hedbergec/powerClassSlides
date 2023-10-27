
sim <- function(delta, n) {
  y.0 <- rnorm(n, mean = 0) #cases in control
  y.1 <- rnorm(n, mean = delta) #cases in treatment
  test <- t.test(y.1,y.0, var.equal = TRUE) #the test
  return(t = test$statistic) #return the result
}

errorPlot <- function(delta, n, alpha = .05, legend = TRUE) {
  df <- 2*n - 2
  ct <- qt(alpha/2, df, lower.tail = FALSE)
  results <- replicate(5000, sim(delta,n))
  ct <- qt(alpha/2, df, lower.tail = FALSE)
  if (delta == 0) {
    histmax <- max(table(round(results , digits = 1)))
    histmin <- min(results)
    ct <- round(qt(alpha/2, df, lower.tail = FALSE), digits = 1)
    sim.power <- 100*
      round((length(results[results >= ct])+
               length(results[results <= -ct]))/
              length(results), digits = 3)
    hist(round(results, digits = 1),
         main = NA,
         xlab = NA,
         col = "white" ,
         breaks = seq(from = min(round(results, digits = 1)),
                      to = max(round(results, digits = 1)), by = .1),
         ylim = c(0,300),xlim = c(-5,5))
    par(new=TRUE)
    hist(round(results[which(round(results, digits = 1) > round(ct, digits = 1) |
                               round(results, digits = 1) < round(-ct, digits = 1)) ],
               digits = 1),
         main = NA,
         xlab = NA,
         col = "black" ,
         breaks = seq(from = min(round(results, digits = 1)),
                      to = max(round(results, digits = 1)),
                      by = .1), 
         ylim = c(0,300),xlim = c(-5,5))
    abline(v = round(c(-ct-.1, ct), digits = 1))
    if (legend == TRUE) {
      legend(-5, 300,
             c(paste("Correct result:"),paste(as.character(100-sim.power),"percent"),
               expression(paste("Type I Error (",alpha,"):")), paste(as.character(sim.power),"percent")),
             fill = c("white", NA, "black", NA), border = c("black", NA, "black", NA), cex =1,box.col = NA)
    }
    
  }
  if (delta != 0)  {
    histmax <- max(table(round(results , digits = 1)))
    histmin <- min(results)
    sim.power <- 100*
      round((length(results[results >= ct])+
               length(results[results <= -ct]))/
              length(results), digits = 3)

    hist(round(results, digits = 1),
         main = NA,
         xlab = NA,
         col = "white" ,
         breaks = seq(from = min(round(results, digits = 1)),
                      to = max(round(results, digits = 1)), by = .1),
         ylim = c(0,300),
         xlim = c(-5,5))

    par(new=TRUE)
    hist(round(results[which(round(results, digits = 1) <= round(ct, digits = 1) & 
                               round(results, digits = 1) >= round(-ct, digits = 1))],
               digits = 1),
         main = NA,
         xlab = paste("Test result, critical value = ",
                      as.character(round(ct, digits = 2))),
         col = "gray" ,
         breaks = seq(from = min(round(results, digits = 1)),
                      to = max(round(results, digits = 1)),
                      by = .1), 
         ylim = c(0,300),
         xlim = c(-5,5))
    abline(v = round(c(-ct-.1, ct), digits = 1))
    if (legend == TRUE) {
      legend(-5, 300,
             c(paste("Correct result:"),paste(as.character(sim.power),"percent"),
               expression(paste("Type II Error (",beta,"):")), paste(as.character(100-sim.power),"percent")),
             fill = c("white", NA, "gray", NA) , border = c("black", NA, "black", NA),  cex = 1, box.col = NA)
    }
    
  }
}

cdfPlot <- function(z) {
  x <- seq(-5, 5, length = 100)
  cden <- dnorm(x,mean = 0)
  z.x <- c(-5, seq(-5, z, 0.01), z)
  z.y <- c(0,dnorm(seq(-5, z, 0.01)),0)
  plot(x, cden, type = "l", ylab = "Density", xlab = "z")
  polygon(z.x, z.y, col = "black")
  lines(x, cden)
  title(main = bquote(paste("z = ",.(z),", ",Phi,"(z)=",.(pnorm(z)))), cex.main = 2)
}

replicationPlot <- function(alpha, beta, proptrue) {
  R <- proptrue/(1-proptrue)
  results <- c((1-beta)*R/(R+1),beta*R/(R+1),(1-alpha)/(R+1),alpha/(R+1))
  cat <- c(paste("Truly significant (published)"),
           paste("Type II Error (not published)"),
           paste("Truly false (not published)"),
           paste("Type I Error (published)"))
  pie(results, 
      labels = cat,
      cex = 1.5,
      col=c("white","gray","white","black"),
      radius = 1) 
  title(main = paste("Percent of significant findings that are true =", 
                      as.character(round(100*(results[1]/(results[1]+results[4]))))), cex.main = 1.75)
}


exactPowerPlot <- function(type, es, n, m = 0, ryx = 0, icc = 0, 
                           r2u = 0, r2c = 0, r2t = 0, upsilon = 0, 
                           cov = "No", alpha = .05,
                           marklambda = FALSE) {
  if (type == "srs") {
    ncp <- es*sqrt(n/2)
    if (cov == "Yes") {
      ncp <- es*sqrt(n/2)*sqrt(1/(1-ryx^2))*sqrt((2*n-3)/(2*n-2))
    }
  }
  else if (type == "crt") {
    ncp <- es*sqrt(n*m/2)*sqrt(1/(1+(n-1)*icc))
    if (cov == "Yes") {
      ncp <- es*sqrt(n*m/2)*sqrt(1/(1+(n-1)*icc-(r2u+(n*r2c-r2u)*icc)))
    }
  }
  else if (type == "msrt") {
    ncp <- es*sqrt(n*m/2)*sqrt(1/(1+(upsilon*n/2-1)*icc))
    if (cov == "Yes") {
      ncp <- es*sqrt(n*m/2)*sqrt(1/(1+(upsilon*n/2-1)*icc-(r2u+(upsilon*n/2*r2t-r2u)*icc)))
    }
  }
  
  if (type == "srs") {
    df <- n*2-2
  }
  else if (type == "crt") {
    df <- m*2-2
  }
  else if (type == "msrt") {
    df <- m-1
  }
  if (cov == "Yes") {
    df <- df-1
  }

  t <- seq(-5, 8, length = 500)
  ncden <- dt(t, df , ncp)
  cden <- dt(t, df)
  critR <- qt(alpha/2, df, lower.tail = FALSE)
  critL <- qt(alpha/2, df, lower.tail = TRUE)
  
  alphaR.x <- c(critR, seq(critR, 5, 0.01), 5)
  alphaR.y <- c(0, dt(seq(critR, 5, 0.01), df), 0)
  
  alphaL.x <- c(critL, seq(-5, critL, 0.01), critL)
  alphaL.y <- c(0, dt(seq(-5, critL, 0.01), df), 0)
  
  beta.x <- c(critL, seq(critL, critR, 0.01), critR)
  beta.y <- c(0, dt(seq(critL, critR, 0.01), df, ncp), 0)
  
  plot(t, cden, type = "l", ylab = "Density", xlab = "t", main = NA)
  polygon(alphaR.x, alphaR.y, col = "black")
  polygon(alphaL.x, alphaL.y, col = "black")
  lines(t, cden)
  
  polygon(beta.x, beta.y, col = "gray", lty = 0)
  lines(t, ncden, lty = 2)
  lines(t, cden)
  if (marklambda == TRUE) {
    arrows(ncp, dt(ncp, df, ncp)-.15, ncp, 0, col = "blue")
    text(ncp, dt(ncp, df, ncp)-.1, expression(lambda), cex = 2)
  }
  beta <- pt(critR,df,ncp) - pt(-critR,df,ncp)
  power <- 1-beta
  legend(-5.5, 0.35, 
         c(expression(paste("Type I Error (",
                                     alpha,")")), 
           expression(paste("Type II Error (",beta,")")),
           bquote(paste("NCP (",lambda,") = ", .(round(ncp, digits = 2)))),
           paste("Power = ", as.character(round(power, digits = 2)))), 
         fill = c("black", "gray", NA, NA), 
         border = "white",
         box.col = NA, cex = 1)
  legend(-5.5, 0.15, 
         c("Null Dist.","Alt. Dist."), 
         lty = c(1,2), 
         box.col = NA, cex = 1)
}

powervalues <- function(n,delta,alpha) {
  df <- 2*n-2
  ct <- qt(alpha/2,df, lower.tail = FALSE)
  ncp <- delta*sqrt(n/2)
  beta <- pt(ct,df,ncp) - pt(-ct,df,ncp)
  power <- 1-beta
  return(power)
}

mdesvalues <- function(n,power,alpha) {
  df <- 2*n-2
  M_t <- qt(power,2*n-2)-qt(alpha/2,2*n-2)
  mdes <- sqrt((2*M_t^2)/n)
  return(mdes)
}


# errorAppUI <- function(id) {
#   ns <- NS(id)
#   inputPanel(
#     sliderInput("n", label = "sample size in each group",
#                 min = 2, max = 50, value = 30, step = 2),
#     sliderInput("delta", label = "effect size",
#                 min = 0, max = 1, value = 0, step = 0.1)
#   )
#   plotOutput(ns("simPlot"))
# }
# 
# errorApp <- function(input, output, session) {
#   output$simPlot <- renderPlot ({
#     errorPlot(delta = input$delta, n = input$n)
#   })
# }
# 
# srsCurveAppUI <- function(id) {
#   ns <- NS(id)
#   inputPanel(
#     sliderInput("n2", label = "Number per group",
#                 min = 2, max = 50, step = 1, value = 25),
#     
#     sliderInput("es", label = "Effect size",
#                 min = 0.1, max = 1, value = .3, step = 0.1)
#     )
#   plotOutput(ns("srsCurvePlot"))
# }
# 
# srsCurveApp <- function(input, output, session) {
#   output$srsCurvePlot <- renderPlot({
#     exactPowerPlot(type = "srs", es = input$es, n = input$n2)
#   })
# }
# 
# replicationPlotAppUI <- function(id) {
#   ns <- NS(id)
#   inputPanel(
#     sliderInput("a", label = "alpha",
#                 min = .001, max = .1, step = .001, value = .05),
#     
#     sliderInput("b", label = "power",
#                 min = 0.2, max = .9, value = .5, step = 0.1),
#     sliderInput("p", label = "Proportion hypotheses tested that are true",
#                 min = 0.025, max = .8, value = .2, step = 0.025)
#   )
#   plotOutput(ns("replicationPlot"))
# }
# 
# replicationPlotApp <- function(input, output, session) {
#   output$replicationPie <- renderPlot({
#     replicationPlot(input$a, 1-input$b, input$p) 
#   })
# }

# computePower <- reactive({
#   req(type)
#   req(cov)
#   req(es)
#   if (type == "srs") {
#     whatisX <- "Units per group"
#     x <- seq(2,50, by=2)
#     N <- 2*x
#     df <- 2*x-2
#     ncp <- es*sqrt(x/2)
#     if (cov == "Yes") {
#       df <- df -1
#       ncp <- es*sqrt(x/2)*sqrt(1/(1-ryx^2))*sqrt((2*x-3)/(2*x-2))
#     } 
#   }
#   else if (type == "crt") {
#     whatisX <- "Clusters per group"
#     x <- seq(2,25)
#     N <- 2*n*x
#     df <- 2*x-2
#     ncp <- es*sqrt(n*x/2)*sqrt(1/(1+(n-1)*icc))
#     if (cov == "Yes") {
#       df <- df - 1
#       ncp <- es*sqrt(n*x/2)*sqrt(1/(1+(n-1)*icc-(r2u+(n*r2c-r2u)*icc)))
#     }
#   }
#   else if (type == "msrt") {
#     whatisX <- "Total clusters"
#     x <- seq(2,25)
#     N <- 2*n*x
#     df <- x-1
#     ncp <- es*sqrt(n*x/2)*sqrt(1/(1+(upsilon*n/2-1)*icc))
#     if (cov == "Yes") {
#       df <- df - 1
#       ncp <- es*sqrt(n*x/2)*sqrt(1/(1+(upsilon*n/2-1)*icc-(r2u+(upsilon*n/2*r2t-r2u)*icc)))
#     }
#   }
#   power <- 1-(pt(qt(alpha/2,df, lower.tail = FALSE),df,ncp) - 
#                 pt(-qt(alpha/2,df, lower.tail = FALSE),df,ncp))
#   return(list(power = power, x = x, df = df, N = N, ncp = ncp, alpha = alpha, type = type, whatisX = whatisX))
# })
# 
# computeMDES <- reactive({
#   req(type)
#   req(cov)
#   req(es)
#   req(power)
#   if (type == "srs") {
#     whatisX <- "Units per group"
#     x <- seq(2,50)
#     N <- 2*x
#     df <- 2*x-2
#     M <- qt(power,df)-qt(alpha/2,df)
#     D <- 1
#     sample <- x
#     if (cov == "Yes") {
#       df <- df - 1
#       M <- qt(power,df)-qt(alpha/2,df)
#       D <- sqrt(1-ryx^2)
#     } 
#   }
#   else if (type == "crt") {
#     whatisX <- "Clusters per group"
#     x <- seq(2,25)
#     N <- 2*n*x
#     df <- 2*x-2
#     M <- qt(power,df)-qt(alpha/2,df)
#     sample <- x*n
#     D <- sqrt(1+(n-1)*icc)
#     
#     if (cov == "Yes") {
#       df <- df - 1
#       M <- qt(power,df)-qt(alpha/2,df)
#       D <- sqrt(1+(n-1)*icc-(r2u+(n*r2c-r2u)*icc))
#     }
#   }
#   else if (type == "msrt") {
#     whatisX <- "Total clusters"
#     x <- seq(2,25)
#     N <- 2*n*x
#     df <- x-1
#     M <- qt(power,df)-qt(alpha/2,df)
#     sample <- x*n
#     D <- sqrt(1+(upsilon*n/2-1)*icc)
#     if (cov == "Yes") {
#       df <- df - 1
#       M <- qt(power,df)-qt(alpha/2,df)
#       D <- sqrt(1+(upsilon*n/2-1)*icc-(r2u+(upsilon*n/2*r2t-r2u)*icc))
#     }
#   }
#   mdes <- M*sqrt(2/sample)*D
#   return(list(MDES = mdes, x = x, df = df, N = N, M = M, power = power, alpha = alpha, type = type, whatisX = whatisX))
# })
# 
# plotMake <- reactive({
#   req(out)
#   req(type)
#   if (out == "ep") {
#     req(length(t) == length(cden))
#     
#     
#   }
#   else if (out == "yPower" & length(computePower()$x) == length(computePower()$power)) {
#     req(length(x) == length(power))
#     return (
#       plot_ly(
#         x = ~computePower()$x, 
#         y = ~computePower()$power,
#         type = 'scatter', 
#         mode = 'lines',
#         hoverinfo = 'text',
#         text = paste("Power=",as.character(round(computePower()$power, digits = 3)), ",",
#                      computePower()$whatisX,"=", computePower()$x)) %>% 
#         add_markers(marker = list(size = 10, color='rgba(80, 80, 80, .5)')) %>%
#         layout (showlegend = FALSE,
#                 xaxis = list(title = computePower()$whatisX),
#                 yaxis = list(title = "Power", range = c(0,1))
#         )
#     )
#   }
#   else if (out == "yMDES" & length(computeMDES()$x) == length(computeMDES()$MDES)) {
#     req(length(x) == length(power)) 
#     return (
#       plot_ly(
#         x = ~computeMDES()$x, 
#         y = ~computeMDES()$MDES,
#         type = 'scatter', 
#         mode = 'lines',
#         hoverinfo = 'text',
#         text = paste("MDES=",as.character(round(computeMDES()$MDES, digits = 3)), ",",
#                      computeMDES()$whatisX,"=", computeMDES()$x)) %>% 
#         add_markers(marker = list(size = 10, color='rgba(80, 80, 80, .5)')) %>%
#         layout (showlegend = FALSE,
#                 xaxis = list(title = computeMDES()$whatisX),
#                 yaxis = list(title = "MDES", range = c(0,3))
#         )
#     )
#   }
#   else {
#     return(
#       plotly_empty(x = 0 , y = 0, text = "loading plot...") %>% add_text(textfont = list(size = 40), textposition = "center")
#     )
#   }
# })
# 
# output$thePlot <- renderPlotly({
#   plotMake()
# })
# 
# output$showThePlot <- renderUI({
#   plotlyOutput("thePlot", height = 400)
# })