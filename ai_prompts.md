# AI Interaction

I used Claude throughout this project to clarify R syntax, learn about the distributions, suggest measures of normality, and help create an interactive sampling distribution explorer. The decisions about what "approximately normal" means and which sample sizes to test were mine.

## Examples

**1. Project setup and structure.**
I asked how an R project is organized and how to push files to GitHub. The answer mixed up scripts and analyses, so I asked again by quoting the assignment ("each analysis should contain R code, results, plots..."). I ended up with shared code in `R/`, one Quarto file per population, a `results/` folder, and a summary file.

**2. Shapiro test.**
Claude suggested a Shapiro–Wilk test and wrote the code for it. I did not understand it, so I asked where it came from and whether I needed it, and told Claude to remove it. I used only the metrics (simulated mean, SD, skewness, kurtosis, theoretical SE) that I understood.

**3. The `flag_good` function.**
I did not understand where `flag_good` came from, so I asked Claude to go through the code line by line. Once I understood that it checks my cutoffs, I decided what the cutoffs should be myself.

**4. Skewness and kurtosis cutoffs.**
I asked what the common standards for skewness and excess kurtosis are so I could be stricter. My first section 6 only stated numbers, so I asked for better defense. I now use |skewness| ≤ 0.25 and |excess kurtosis| ≤ 0.5, and each is compared with the Normal control case.

**5. Formatting problems.**
I asked for help debugging formatting, for example the boxed answer for n. The local preview showed raw `$$\boxed{...}$$` text. I learned it is only drawn on GitHub, so I checked it there.

## How I checked the AI's work

I rendered every file and compared the numbers with theory instead of trusting the text. For example, the die's simulated excess kurtosis at n = 1 was −1.268 against a theoretical −1.269. I re-rendered earlier files after each change to shared functions to confirm the results did not change, and I read each interpretation against its table. I went through and understood the code before writing the summary.