# Module 2: Git and GitHub for Beginners

Have you ever saved files like this?
*   `analysis.R`
*   `analysis_final.R`
*   `analysis_final_v2.R`
*   `analysis_really_final.R`

This is a messy way to track changes. **Git** is a tool that solves this problem. It works like a time machine for your project, keeping a neat history of every change you make to your files.

In this chapter, we will learn how Git tracks files and how **GitHub** lets you share those files online.

---

## 1. What is Git & How Does it Work?

Git divides your project into three areas:

```
  [ Working Directory ]    --->    [ Staging Area ]    --->    [ Local Repository ]
  (Where you edit files)           (Where you stage)            (Where Git saves permanently)
```

1.  **Working Directory**: The folder on your computer where you edit files.
2.  **Staging Area (Index)**: A "loading dock" where you select which changes you want to save.
3.  **Local Repository**: The database where Git permanently saves your changes as "checkpoints" called **commits**.

---

## 2. Basic Git Commands

Here are the basic commands you will run in your terminal.

### A. Start Tracking a Folder (`git init`)
To tell Git to start watching a folder, go to the folder in your terminal and type:
```fish
git init
```

### B. Checking Status (`git status`)
To see which files have changed and what Git is currently watching, type:
```fish
git status
```

### C. Staging a File (`git add`)
To add a file to the "staging area" (loading dock), use `git add`:
```fish
# Stage a specific file
git add script.R

# Stage all files in the current folder
git add .
```

### D. Saving Your Changes (`git commit`)
To save a snapshot checkpoint of your staged files, run `git commit` with a brief message explaining what you changed:
```fish
git commit -m "Create initial script for loading sample data"
```

---

## 3. Remote Repositories & GitHub

**GitHub** is a website where you can upload your local Git history. This lets you collaborate with others or keep a backup of your code.

*   `git push`: Sends your local saves up to GitHub.
*   `git pull`: Downloads files and updates from GitHub to your computer.

```fish
# 1. Connect your local folder to a GitHub repository link
git remote add origin https://github.com/yourusername/my_project.git

# 2. Upload your files to GitHub
git push -u origin main

# 3. Download updates from GitHub
git pull origin main
```

### A. Authentication Security: PAT & SSH Keys
To interact with GitHub securely, passwords are not allowed. You must use one of these two methods:
1.  **Personal Access Tokens (PAT)**: Generated in GitHub Developer Settings. Paste this token when Git prompts you for a password in the terminal.
2.  **SSH Key Pairs**: Public/private key pairs that securely authenticate your local computer to GitHub without typing passwords.
    ```fish
    # Generate a secure SSH key pair
    ssh-keygen -t ed25519 -C "your_email@example.com"
    
    # Add your private key to the SSH agent
    eval (ssh-agent -c)
    ssh-add ~/.ssh/id_ed25519
    ```
    Copy your public key (`cat ~/.ssh/id_ed25519.pub`) and save it in your GitHub settings under "SSH and GPG keys".

---

## 4. Branching & Collaboration Strategy

> **The Problem:** You and a colleague are working on the same script. You both make edits at the same time on the same line, overwriting each other's changes and breaking the main pipeline.
>
> **The Solution:** **Branching** lets you work in separate copy universes of the code, merging them only after they are tested and verified.

When working on a team project, the `main` (or `master`) branch is kept for stable, clean code. All new features or bioinformatic scripts are developed in separate, isolated **branches** to avoid breaking the core project until they are tested.

```
                       +-- [ Feature Branch ] -- (Develop & Test) --+
                       |                                            |
  === [ Main Branch ] =+============================================+=> (Merge)
```

*   **HEAD Pointer**: The invisible pointer that tells Git which commit/branch your working directory is currently viewing.
*   `git branch [name]`: Creates a new branch with a specific name.
*   `git checkout [name]`: Switches your working directory to the target branch.
*   `git merge [name]`: Combines changes from the target branch into your active branch.

### Typical Branching Workflow
```fish
# 1. Create and switch to a new branch for UMAP plotting
git branch feat_umap
git checkout feat_umap

# 2. Edit, stage, and commit changes on the branch
git add umap_plot.R
git commit -m "Create UMAP plotting script"

# 3. Switch back to the main branch
git checkout main

# 4. Merge UMAP script into main
git merge feat_umap
```

### Branch Cleanup & Deletion
Once feature branches are merged and verified, clean them up:
```fish
# Delete branch locally
git branch -d feat_umap

# Delete branch on GitHub remote repository
git push origin -d feat_umap
```

### C. Resolving Merge Conflicts

> **The Problem:** You try to merge a feature branch, but the operation crashes with: `CONFLICT (content): Merge conflict in analysis.R. Automatic merge failed`. This happens because you and a colleague edited the exact same line of code in different ways on different branches, and Git doesn't know which version is correct.
>
> **The Solution:** **Manual Conflict Resolution**. Git marks the conflicted lines inside the file using conflict markers. You must edit the file to select the correct code, remove the markers, and commit the resolved code.

#### How to Read Conflict Markers
When a conflict occurs, open the affected file. You will see markers like this:

```R
<<<<<<< HEAD
# Your code on the active branch (e.g., main)
tpm_normalized <- log2(counts + 0.1)
=======
# Code from the branch you are merging (e.g., feat_umap)
tpm_normalized <- log10(counts + 1)
>>>>>>> feat_umap
```

*   `<<<<<<< HEAD`: Marks the start of your active branch's changes.
*   `=======`: Separates the two conflicting versions.
*   `>>>>>>> feat_umap`: Marks the end of the incoming branch's changes.

#### Resolution Steps
1.  **Open the file** and discuss with your colleague to decide which version of the code is correct.
2.  **Edit the file** to delete the conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`) and keep only the correct line of code.
3.  **Stage and Commit** the resolved file:
    ```fish
    git add analysis.R
    git commit -m "Resolve merge conflict in analysis.R, using log2 normalization"
    ```

---

## 5. Keeping Your Project Safe: The `.gitignore` File

> **The Problem:** You try to push your project folder to GitHub, but it gets blocked with: `Error: File data.bam is 250 MB; this exceeds GitHub's file size limit of 100 MB`. Now your Git history is bloated with heavy files you cannot easily upload.
>
> **The Solution:** The **`.gitignore` file** acts as a filter, telling Git to ignore large data files and only track your code scripts.

When doing bioinformatics, your data files (like large gene spreadsheets or images) are huge. 
*   GitHub does not allow files larger than 100MB.
*   More importantly, you should **never** track raw data files in Git. You only want to track your code.

To tell Git to ignore these heavy files, create a text file named exactly `.gitignore` in your project folder, and list the file patterns you want to ignore.

### Example `.gitignore` File Content
```gitignore
# Ignore large Excel and CSV spreadsheets
*.csv
*.xlsx
*.tsv

# Ignore large sequencing and alignment data formats
*.fastq.gz
*.bam
*.sam

# Ignore hidden folders
.venv/
.DS_Store
```

---

## 🧑‍💻 Hands-On Practice Examples & Exercises

### Scenario: Creating Your First Repository
Let's practice the basic Git commands in your terminal (using your Fish shell):

1.  Create a test folder on your computer and navigate into it:
    ```fish
    mkdir -p ~/git_practice
    cd ~/git_practice
    ```
2.  Initialize the folder as a Git repository.
3.  Create a blank file named `data_cleaner.R`:
    ```fish
    echo "# R code to clean data" > data_cleaner.R
    ```
4.  Run the command to check the status of your folder. Note that `data_cleaner.R` is marked as an "untracked file" (in red).
5.  Stage `data_cleaner.R` using the `git add` command.
6.  Check the status again. It should now be in green ("changes to be committed").
7.  Commit the file with a clear, descriptive message (like `"Add initial cleaning script"`).

---

### Solutions & Expected Terminal Outputs

<details>
<summary>Click to view expected outputs</summary>

#### Command sequence:
```fish
# Initialize repository
git init
# Output: Initialized empty Git repository in /home/user/git_practice/.git/

# Check status
git status
# Output shows data_cleaner.R under "Untracked files"

# Stage the file
git add data_cleaner.R

# Verify status
git status
# Output shows "Changes to be committed: new file: data_cleaner.R"

# Commit changes
git commit -m "Add initial cleaning script"
# Output: [main (root-commit) ...] Add initial cleaning script
```
</details>
