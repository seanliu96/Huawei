echo "##Author: Sean Liu" > README.md
date | xargs -I{} echo "##TIME:" {} >> README.md
echo "###然后我安静的发现
###两个人已经没有任何语言
###曾经你纯真的永远
###让我不顾一切开始怀念" >> README.md

git add `git ls-files -o`
git add `git ls-files -m`
git rm `git ls-files -d`
date | xargs -I{} git commit -m {}
git push origin master
