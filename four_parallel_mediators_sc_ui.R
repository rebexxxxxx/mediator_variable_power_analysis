#three parallel mediators sc UI
list(
  withTags(
    table(style = "width: 350px;", id = "STDpath_table",
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a1 low")),  
            td(textInput(inputId = "a1_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a1 high")),  
            td(textInput(inputId = "a1_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a2 low")),  
            td(textInput(inputId = "a2_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a2 high")),  
            td(textInput(inputId = "a2_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a3 low")),  
            td(textInput(inputId = "a3_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a3 high")),  
            td(textInput(inputId = "a3_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a4 low")),  
            td(textInput(inputId = "a4_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a4 high")),  
            td(textInput(inputId = "a4_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b1 low")),  
            td(textInput(inputId = "b1_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b1 high")),  
            td(textInput(inputId = "b1_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b2 low")),  
            td(textInput(inputId = "b2_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b2 high")),  
            td(textInput(inputId = "b2_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b3 low")),  
            td(textInput(inputId = "b3_low", label = NULL, value = "0.14")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b3 high")),  
            td(textInput(inputId = "b3_high", label = NULL, value = "0.39"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b4 low")),  
            td(textInput(inputId = "b4_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b4 high")),  
            td(textInput(inputId = "b4_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("c low")),  
            td(textInput(inputId = "c_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("c high")),  
            td(textInput(inputId = "c_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm1m2 low")),  
            td(textInput(inputId = "rm1m2_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm1m2 high")),  
            td(textInput(inputId = "rm1m2_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm1m3 low")),  
            td(textInput(inputId = "rm1m3_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm1m3 high")),  
            td(textInput(inputId = "rm1m3_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm1m4 low")),  
            td(textInput(inputId = "rm1m4_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm1m4 high")),  
            td(textInput(inputId = "rm1m4_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm2m3 low")),  
            td(textInput(inputId = "rm2m3_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm2m3 high")),  
            td(textInput(inputId = "rm2m3_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm2m4 low")),  
            td(textInput(inputId = "rm2m4_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm2m4 high")),  
            td(textInput(inputId = "rm2m4_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm3m4 low")),  
            td(textInput(inputId = "rm3m4_low", label = NULL, value = "0.04")),
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("rm3m4 high")),  
            td(textInput(inputId = "rm3m4_high", label = NULL, value = "0.19"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. dev: X")),  
            td(textInput(inputId = "SDX", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. dev: M1")),  
            td(textInput(inputId = "SDM1", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. dev: M2")),  
            td(textInput(inputId = "SDM2", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. dev: M3")),  
            td(textInput(inputId = "SDM3", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. dev: M4")),  
            td(textInput(inputId = "SDM3", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. dev: Y")),  
            td(textInput(inputId = "SDY", label = NULL, value = "1.00"))
          )
    )  
  )
)