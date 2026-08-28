# Module 3: Good Coding Practices

Imagine opening a script you wrote six months ago, only to see variables named `x`, `temp`, `data2`, and lines of code piled together without a single comment. It doesn't run, and you have absolutely no idea what it is supposed to do. Now imagine a colleague trying to use that script to replicate your publication results.

This is the cost of "spaghetti code." **Good Coding Practices** are rules that keep your code clean, readable, and reproducible. In this chapter, we will learn how to write neat, professional scripts that you and your lab members can easily understand years from now.

---

## 1. Naming Things Clearly

When you write code, you create variables, functions, and files. Choosing good names is one of the most important habits to build.

### A. Case Styles
Choose one style and stick to it. In the R coding community, **snake_case** is the most popular:
*   `snake_case`: `my_data_table`, `cell_counts` (Recommended)
*   `camelCase`: `myDataTable`, `cellCounts`

### B. General Naming Tips
1.  **Variables**: Use clear nouns that describe what is stored.
    *   *Bad*: `x`, `d`, `data`, `v1`
    *   *Good*: `raw_expression_counts`, `sample_metadata`, `filtered_spots`
2.  **Functions**: Use action verbs that state what the function does.
    *   *Bad*: `my_func`, `calculation`, `normalization`
    *   *Good*: `calculate_average`, `filter_failed_samples`
3.  **Files**: Prefix files with numbers if they run in a specific order. Avoid spaces, capital letters, and special characters.
    *   *Bad*: `Analysis Script final version 2.R`
    *   *Good*: `01_read_data.R`, `02_clean_data.R`, `03_plot_results.R`

---

## 2. What is a Function & Why Use Them? (Modularity)

Imagine you have to convert temperatures from Fahrenheit to Celsius multiple times in your analysis. Instead of copying and pasting the math formula over and over, you write a **function** once and reuse it.

This is called the **DRY (Don't Repeat Yourself)** principle.

### Writing a Simple Function
Here is how you write a custom function in R:

```R
# Define a function named fahrenheit_to_celsius
fahrenheit_to_celsius <- function(temp_f) {
  # Perform the math conversion
  temp_c <- (temp_f - 32) * 5 / 9
  
  # Return the final answer
  return(temp_c)
}

# Now we can reuse it easily!
print(fahrenheit_to_celsius(32))  # Output: 0
print(fahrenheit_to_celsius(212)) # Output: 100
```

---

## 3. Documenting Your Code

Documentation means writing clear explanations of what your code is doing. In R, we use `#` comments to document functions.

### Simple Function Documentation
A well-documented function should explain:
1. What the function does.
2. What input values it expects.
3. What output value it returns.

```R
# This function calculates the average count of cells in a data frame column.
#
# Arguments:
#   dataset: a data frame containing sample records
#   column_name: a character string representing the column to average
#
# Returns:
#   A single numeric value representing the average cell count.
calculate_average_cells <- function(dataset, column_name) {
  values <- dataset[, column_name]
  average <- mean(values, na.rm = TRUE)
  return(average)
}
```

---

## 4. Defensive Programming (Preventing Errors)

Defensive programming means checking that your inputs are valid *before* performing calculations. This helps prevent mysterious crashes and tells the user exactly what went wrong.

### A. Using `stop()`
The `stop()` function halts code execution and prints an error message.

```R
calculate_percentage <- function(count, total) {
  # Stop the code if the total is zero (to prevent dividing by zero)
  if (total == 0) {
    stop("Error: Total cannot be zero!")
  }
  
  percentage <- (count / total) * 100
  return(percentage)
}

# This will run fine:
calculate_percentage(5, 20) # Output: 25

# This will trigger our error message and stop running:
calculate_percentage(5, 0)
# Output: Error: Total cannot be zero!
```

---

## 🧑‍💻 Hands-On Practice Examples & Exercises

### Refactoring Challenge: Cleaning Up Messy Code
Below is a messy script written by a researcher. It has poor variable names, hardcoded calculations, and no comments.

```R
# messy.R
d <- data.frame(val = c(12, 18, 5, 22))
v1 <- d$val[1] * 2
v2 <- d$val[2] * 2
v3 <- d$val[3] * 2
v4 <- d$val[4] * 2
res <- c(v1, v2, v3, v4)
print(res)
```

### Tasks:
1.  **Refactor this script** to make it clean and modular.
2.  Write a reusable function called `double_values` that:
    *   Accepts a numeric vector as input.
    *   Checks if the input is empty or not numeric using `is.numeric()`, and stops with a friendly error if it's incorrect.
    *   Multiplies the vector elements by 2.
    *   Returns the doubled vector.
3.  Add clear comments explaining what the function does.

---

### Solutions

<details>
<summary>Click to view solutions</summary>

#### Refactored Script (`clean_script.R`)

```R
# This function takes a list of numbers and doubles each value.
#
# Arguments:
#   numbers: a vector of numeric values
#
# Returns:
#   A numeric vector where each element is multiplied by 2.
double_values <- function(numbers) {
  # Check if input is actually numbers
  if (!is.numeric(numbers)) {
    stop("Error: Input must be a numeric vector.")
  }
  
  # Multiply all numbers by 2
  doubled_result <- numbers * 2
  
  return(doubled_result)
}

# --- Execution ---
# Define sample data frame with clear names
expression_data <- data.frame(counts = c(12, 18, 5, 22))

# Use our function to double the counts column
doubled_counts <- double_values(expression_data$counts)

print(doubled_counts)
# Output: 24 36 10 44
```
</details>
