# Some command lines for cloning, creating branch and pushing

# https://dont-be-afraid-to-commit.readthedocs.io/en/latest/git/commandlinegit.html#some-basic-git-operations # Explanation of some of following command can be found in this website.

$ cd /idia/users/aycha

$ git clone https://github.com/idia-astro/gps-mightee.git # clone gps-mightee to local, get a directory of gps-mightee 

$ cd /idia/users/aycha/gps-mightee

$ git branch -a # show all the branch

$ git checkout -b Sambatra_branch master # create a new branch Sambatra_branch, checkout to master
                                         # There will be an asterisk (*) next to the branch that you are currently on.

$ git push -u origin Sambatra_branch     # push the new created branch into github.

$ git init                               # To initialize in the custom branch

$ git add *                              # To add files in the custom branch

# Before committing, you need to specify who you are:

$ git config --global user.email "aychasam@gmail.com"

$ git config --global user.name "aychasam"

$ git commit -m "Files pushed into Sambatra_branch" # To commit the changes made in the custom branch into a message

$ git push                               # To make changes in your GitHub repo

# If a file was modified

# Modify the file

$ git add <file>

$ git commit -m "git_steps.md modified"                                      # always commit before pushing

$ git push -u origin Sambatra_branch
