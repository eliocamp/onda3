# Add an id attribute based on the label
old_label <- knitr::opts_hooks$get("label")
knitr::opts_hooks$set(label = function(options) {
   if (!is.null(old_label)) {
      options <- old_label(options)
   }
   if (is.null(options[["out.extra"]])) {
      options[["out.extra"]] <- paste0("id='", options[["label"]], "'")
   }
   options
})


old_plot <- knitr::knit_hooks$get("plot")
knitr::knit_hooks$set(plot = function(x, options) {
   # Only show the first plot
   file <- xfun::sans_ext(basename(x))
   number <- as.numeric(rev(strsplit(file, "-", TRUE)[[1]])[1])
   if (!file.exists("text.txt")) {
      writeLines(file, "test.txt")   
   }
   
   x <- old_plot(x, options)
   # print(options)
   if (is.null(options[["one_figure"]])) {
      return(x)
   }
   
   if (number > 1) {
      return("")
   }
   
   # Add the plot-switching script to the document below the plot
   id <- options$label
   
   id_safe <- gsub("-", "_", id)
   id_safe <- gsub("\\.", "_", id_safe)
   script <- glue::glue('
<script>
function change_{id_safe}(newVal){{
  document.getElementById("{id}").src = "{knitr::opts_chunk$get("fig.path")}/{id}-" + newVal + ".png";
}}
</script>
')
   
   script_self <- glue::glue(
      '<script>
      function cycle_{id_safe}(){{
         image = document.getElementById("{id}");
         if (!image.hasAttribute("data-current-img")) {{
            current = 1;
         }} else {{
            current = image.getAttribute("data-current-img");
         }}
         
         if (current == 2) {{
            image.src = "{knitr::opts_chunk$get("fig.path")}/{id}-" + 1 + ".png";
            image.setAttribute("data-current-img", 1);
         }} else {{
            image.src = "{knitr::opts_chunk$get("fig.path")}/{id}-" + 2 + ".png";
            image.setAttribute("data-current-img", 2);
         }}
      
   }}
      </script>'
   )
   # writeLines(x, options$label)
   # pre <- glue::glue('<a href ="javascript:cycle_{id_safe}()"><div class="figure">')
   # 
   # x <- gsub('<div class="figure">', pre, x)
   # x <- gsub('</div>', '</div></a>', x)
   # 
   
   paste0(x, "\n", script, "\n")
}
)


# Customise the knitr_print function for lists
# to supress the [1], [2], etc... 
knitr_print.list <- function(x, options, ...) {
   if (inherits(x[[1]], "gg")) { # if the list has ggplot objects
      for (i in seq_along(x)) {
         knitr::knit_print(x[[i]])
      }  
   } else {
      NextMethod("knitr_print")
   }
}
library(knitr)
registerS3method("knit_print", "list", knitr_print.list)

# Helper functions. 

# Add this as a link to create an hyperlink that changes the plot
# to a specific plot
change_plot <- function(id, value) {
   id_safe <- gsub("-", "_", id)
   id_safe <- gsub("\\.", "_", id_safe)
   knitr::asis_output(glue::glue("javascript:change_{id_safe}('{value}')"))
}

# Add a series of buttons to change the plot
change_plot_buttons <- function(id, names, values = seq_along(names)) {
   id_safe <- gsub("-", "_", id)
   id_safe <- gsub("\\.", "_", id_safe)
   
   lapply(seq_along(values), function(i) {
      value <- values[i]
      name <- names[i]
      knitr::asis_output(glue::glue("<button type=\"button\" onclick=\"javascript:change_{id_safe}('{value}')\">{name}</button>"
      )
      )
   }
   )
}

change_plot_links <- function(id, names, values = seq_along(names)) {
   id_safe <- gsub("-", "_", id)
   id_safe <- gsub("\\.", "_", id_safe)
   
   words <- vapply(seq_along(values), function(i) {
      value <- values[i]
      name <- names[i]
      glue::glue("<a href=\"javascript:change_{id_safe}('{value}')\">{name}</a>")
   }, character(1)
   )
   
   words <- paste0(words, collapse = ", ")
   knitr::asis_output(words)
}
