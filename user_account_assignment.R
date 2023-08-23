#####管理员为用户生成User Account(需要在服务器上跑这个代码)#####
##从1-1000000中随机抽取一个数字作为user account
#Basic Settings
##设置工作目录
# setwd("C:\\Users\\LEGION\\Desktop\\MOVICS_RShiny\\workingDirectory") #本地测试目录
# setwd("/shinyapps/MOVICS_RShiny/workingDirectory") #服务器测试目录
setwd("/srv/shiny-server/workingDirectory") #服务器测试目录
workdir <- getwd()
sam <- sample(1:1000000,1)
while(dir.exists(file.path(workdir,sam))){
  sam <- sample(1:1000000,1)
}
print(sam)
#####Username默认使用用户发送申请邮件的邮箱,将Username和User Account写到"usersInfo.txt"中,分别对应"user_name"和"user_account"列,然后将"usersInfo.txt"文件覆盖上传到"/srv/shiny-server/workingDirectory"路径下#####