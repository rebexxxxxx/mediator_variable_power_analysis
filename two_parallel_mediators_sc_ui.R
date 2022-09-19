#two parallel mediator variables std coef ui
list(
  withTags(
    table(style = "width: 150px;", id = "STDpath_table",
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a1")),  
            td(textInput(inputId = "STa1", label = NULL, value = "0.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("a2")),  
            td(textInput(inputId = "STa2", label = NULL, value = "0.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b1")),  
            td(textInput(inputId = "STb1", label = NULL, value = "0.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("b2")),  
            td(textInput(inputId = "STb2", label = NULL, value = "0.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("c'")),  
            td(textInput(inputId = "STcprime", label = NULL, value = "0.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("r")),  
            td(textInput(inputId = "cor32", label = NULL, value = "0.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. Deviation: X")),  
            td(textInput(inputId = "SDX", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. Deviation: M1")),  
            td(textInput(inputId = "SDM1", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. Deviation: M2")),  
            td(textInput(inputId = "SDM2", label = NULL, value = "1.00"))
          ),
          tr(
            td(style="padding-top:0px;padding-left:10px;width:32px;", label("std. Deviation: Y")),  
            td(textInput(inputId = "SDY", label = NULL, value = "1.00"))
          )
    )  
  )
)