# Some command lines for cloning, creating a custom branch and pushing it to the repository
Explanation of some of following command can be found in this website:
https://dont-be-afraid-to-commit.readthedocs.io/en/latest/git/commandlinegit.html#some-basic-git-operations 

## Cloning
`cd /idia/users/aycha`\
Clone gps-mightee to local, get a directory of gps-mightee\
`git clone https://github.com/idia-astro/gps-mightee.git`\
`cd /idia/users/aycha/gps-mightee`\
`git branch -a # show all the branch`\
## New branch on local directory
Create a new branch Sambatra_branch, checkout to master\
`git checkout -b Sambatra_branch master`\
There will be an asterisk (*) next to the branch that you are currently on.\
## New branch appearing on github
Push the new created branch into github:\
`git push -u origin Sambatra_branch`\
To initialize in the custom branch:\
`git init`\
## Adding files on github
To add files in the custom branch:\
`git add *` or `git add file`\                   
## Before committing, you need to specify who you are:
`git config --global user.email "aychasam@gmail.com"`\
`git config --global user.name "aychasam"`\
To commit the changes made in the custom branch into a message (N.B: always commit before push)\
`git commit -m "Files pushed into Sambatra_branch"`\
To make changes in your customized branch in GitHub repository. So, if you want your files and commits to appear in the master branch, before commiting you can exit your branch and go to master with: `git switch master`. But if you want modifications only in your branch, type:\
`git push -u origin Sambatra_branch` 
## If a file was modified and you want to add it to your branch
Modify the file, then:\
`git add fileName`\
`git commit -m "git_steps.md modified"`\# always commit before pushing
`git push -u origin Sambatra_branch`
## Pull and merge a file into the master branch