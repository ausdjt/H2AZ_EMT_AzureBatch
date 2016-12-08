
ui = shinyUI(fluidPage(
  h3("MetricsGraphics Example", style="text-align:center"),
  metricsgraphicsOutput('mjs1')
))

  
  tgfb <- emtSignatureData[["conditionMDCKTGFb"]]$dataTable %>%
    mjs_plot(x=b, y=log10_qval, linked=TRUE, width=350, height=275)
  
  shZ <- emtSignatureData[["conditionMDCKshZ"]]$dataTable %>%
    mjs_plot(x=b, y=log10_qval, linked=TRUE, width=350, height=275)
  
  output$mjs1 <- mjs_grid(tgfb, shZ, ncol=2)
  
  

  tmp <- data.frame(year=seq(1790, 1970, 10), uspop=as.numeric(uspop))
  
  tmp %>% mjs_plot(x=year, y=uspop) %>% mjs_line() %>%mjs_add_marker(1850, "Something Wonderful") %>%  mjs_add_baseline(150, "Something Awful")
  
  
  set.seed(1492)
  stocks <- data.frame(
    time = as.Date('2009-01-01') + 0:9,
    X = rnorm(10, 0, 1),
    Y = rnorm(10, 0, 2),
    Z = rnorm(10, 0, 4))
  
  stocks2 <- data.frame(
    time = as.Date('2009-01-01') + 0:9,
    X = rnorm(10, 0, 1),
    Y = rnorm(10, 0, 2),
    Z = rnorm(10, 0, 4))
  
  stocks %>%
    mjs_plot(x=time, y=X) %>%
    mjs_line() %>%
    mjs_add_line(Y) %>%
    mjs_add_line(Z) %>%
    mjs_axis_x(xax_format="date") %>%
    mjs_add_legend(legend=c("X", "Y", "Z"))
  
  s1 <- stocks %>%
    mjs_plot(x=time, y=X, linked=TRUE, width=350, height=275) %>%
    mjs_line() %>%
    mjs_add_line(Y) %>%
    mjs_add_line(Z) %>%
    mjs_axis_x(xax_format="date") %>%
    mjs_add_legend(legend=c("X", "Y", "Z"))
  
  s2 <- stocks2 %>%
    mjs_plot(x=time, y=X, linked=TRUE, width=350, height=275) %>%
    mjs_line() %>%
    mjs_add_line(Y) %>%
    mjs_add_line(Z) %>%
    mjs_axis_x(xax_format="date") %>%
    mjs_add_legend(legend=c("X", "Y", "Z"))

html_print(  mjs_grid(s1, s2, ncol=2)  )

mtcars %>%
  mjs_plot(x=wt, y=mpg, width=600, height=500) %>%
  mjs_point(least_squares=TRUE) %>%
  mjs_labs(x="Weight of Car", y="Miles per Gallon")

  