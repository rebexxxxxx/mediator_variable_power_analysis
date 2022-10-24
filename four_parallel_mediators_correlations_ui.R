#three parallel mediators corr UI
list(
  withTags(
    table(id = "correlations_table", style = "width:500px;",
          tr(
            th(style="padding-top:5px;padding-left:10px;width:40px;"),
            th(style="padding-top:5px;", label("X Low")),
            th(style="padding-top:5px;",label("M1 Low")),
            th(style="padding-top:5px;",label("M2 Low")),
            th(style="padding-top:5px;",label("M3 Low")),
            th(style="padding-top:5px;",label("M4 Low")),
            th(style="padding-top:5px;",label("Y Low"))
          ),
          tr(
            td(style="padding-left:10px;width:40px;", label("X Low")),
            td(style="text-align:center", "1.00"),
            td(),
            td(),
            td(),
            td(),
            td()
          ),
          tr(
            td(style="padding-left:10px;width:40px;", label("M1 Low")),
            td(textInput(inputId = "cor21_low", label = NULL, value = "0.14")),
            td(style="text-align:center", "1.00"),
            td(),
            td(),
            td(),
            td()
          ),
          tr(
            td(style="padding-left:10px;width:40px;", label("M2 Low")),  
            td(textInput(inputId = "cor31_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor32_low", label = NULL, value = "0.14")),
            td(style="text-align:center", "1.00"),
            td(),
            td(),
            td()
          ),
          tr(
            td(style="padding-left:10px;width:40px;", label("M3 Low")),  
            td(textInput(inputId = "cor41_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor42_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor43_low", label = NULL, value = "0.14")),
            td(style="text-align:center;", "1.00"),
            td(),
            td()
          ),
          tr(
            td(style="padding-left:10px;width:40px;", label("M4 Low")),  
            td(textInput(inputId = "cor51_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor52_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor53_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor54_low", label = NULL, value = "0.14")),
            td(style="text-align:center;", "1.00"),
            td()
          ),
          tr(
            td(style="padding-left:10px;width:40px;", label("Y Low")),  
            td(textInput(inputId = "cor61_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor62_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor63_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor64_low", label = NULL, value = "0.14")),
            td(textInput(inputId = "cor65_low", label = NULL, value = "0.14")),
            td(style="text-align:center;", "1.00")
          ),
          table(id = "correlations_table", style = "width:500px;",
                tr(
                  th(style="padding-top:5px;padding-left:10px;width:40px;"),
                  th(style="padding-top:5px;", label("X High")),
                  th(style="padding-top:5px;",label("M1 High")),
                  th(style="padding-top:5px;",label("M2 High")),
                  th(style="padding-top:5px;",label("M3 High")),
                  th(style="padding-top:5px;",label("M4 High")),
                  th(style="padding-top:5px;",label("Y High"))
                ),
                tr(
                  td(style="padding-left:10px;width:40px;", label("X High")),
                  td(style="text-align:center", "1.00"),
                  td(),
                  td(),
                  td(),
                  td(),
                  td()
                ),
                tr(
                  td(style="padding-left:10px;width:40px;", label("M1 High")),
                  td(textInput(inputId = "cor21_high", label = NULL, value = "0.39")),
                  td(style="text-align:center", "1.00"),
                  td(),
                  td(),
                  td(),
                  td()
                ),
                tr(
                  td(style="padding-left:10px;width:40px;", label("M2 High")),  
                  td(textInput(inputId = "cor31_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor32_high", label = NULL, value = "0.39")),
                  td(style="text-align:center", "1.00"),
                  td(),
                  td(),
                  td()
                ),
                tr(
                  td(style="padding-left:10px;width:40px;", label("M3 High")),  
                  td(textInput(inputId = "cor41_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor42_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor43_high", label = NULL, value = "0.39")),
                  td(style="text-align:center;", "1.00"),
                  td(),
                  td()
                ),
                tr(
                  td(style="padding-left:10px;width:40px;", label("M4 High")),  
                  td(textInput(inputId = "cor51_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor52_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor53_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor54_high", label = NULL, value = "0.39")),
                  td(style="text-align:center;", "1.00"),
                  td()
                ),
                tr(
                  td(style="padding-left:10px;width:40px;", label("M4 High")),  
                  td(textInput(inputId = "cor61_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor62_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor63_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor64_high", label = NULL, value = "0.39")),
                  td(textInput(inputId = "cor65_high", label = NULL, value = "0.39")),
                  td(style="text-align:center;", "1.00")
                ),
          
          tr(
            td(style="padding-left:10px;width:40px;", label("Std. Deviation")),  
            td(textInput(inputId = "SDX", label = NULL, value = "1.00")),
            td(textInput(inputId = "SDM1", label = NULL, value = "1.00")),
            td(textInput(inputId = "SDM2", label = NULL, value = "1.00")),
            td(textInput(inputId = "SDM3", label = NULL, value = "1.00")),
            td(textInput(inputId = "SDM4", label = NULL, value = "1.00")),
            td(textInput(inputId = "SDY", label = NULL, value = "1.00"))
          )
    )  
  )
)
)