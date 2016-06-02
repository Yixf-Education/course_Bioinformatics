%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Title: 天津医科大学（TIJMU）实验课教案
%Version: 2013-09-17
%By: 伊现富(Yixf) 
%Email: yixf1986@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\LoadClass[12pt,a4paper]{article}


\RequirePackage[top=0.6in,bottom=0.4in,left=0.6in,right=0.6in]{geometry}

\RequirePackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Courier New}

\RequirePackage[BoldFont,SlantFont,CJKchecksingle,CJKnumber]{xeCJK}
\setCJKmainfont[BoldFont={Adobe Heiti Std},ItalicFont={Adobe Kaiti Std}]{Adobe Song Std}
\setCJKsansfont{Adobe Heiti Std}
\setCJKmonofont{Adobe Fangsong Std}
\punctstyle{hangmobanjiao}

\defaultfontfeatures{Mapping=tex-text}
\RequirePackage{xunicode}
\RequirePackage{xltxtra}

\XeTeXlinebreaklocale "zh"
\XeTeXlinebreakskip = 0pt plus 1pt minus 0.1pt

\linespread{1}

%\RequirePackage{indentfirst}
\makeatletter
\let\@afterindentfalse\@afterindenttrue
\@afterindenttrue
\makeatother
\setlength{\parindent}{2em}

%插入图片
\usepackage{graphicx}
\graphicspath{{figures/}}
%图文混排
%\usepackage{wrapfig}
%\usepackage{picinpar}
\usepackage{picins}
\usepackage{floatflt}

%获取总页数
\RequirePackage{lastpage}

%调整列表间及列表前后的行距
\RequirePackage{enumitem}
\setlist{parsep=0pt,itemsep=2pt,topsep=2pt,partopsep=0pt}
\setlist[enumerate]{widest=0}

%分栏
\RequirePackage{multicol}

%下划线
\RequirePackage{ulem}

%表格
\RequirePackage{multirow}
\RequirePackage{tabu}
\RequirePackage{booktabs}

%行距
\RequirePackage{setspace}

%常用模块
\RequirePackage{xcolor}
\RequirePackage{hyperref}
%\hypersetup{xetex,bookmarksnumbered=true,bookmarksopen=true,pdfborder=1,breaklinks,colorlinks,linkcolor=cyan,filecolor=black,urlcolor=blue,citecolor=green}
\hypersetup{xetex,bookmarksnumbered=true,bookmarksopen=true,pdfborder=1,breaklinks}
\RequirePackage{amsmath,amssymb}


\newcommand\kch{\textbf{课程名称：}}
\newcommand\kecheng[1]{\def\kchstring{#1}}
\newcommand\shy{\textbf{实验名称：}}
\newcommand\shiyan[1]{\def\shystring{#1}}
\newcommand\jsh{\textbf{教师姓名：}}
\newcommand\jiaoshi[1]{\def\jshstring{#1}}
\newcommand\zhch{\textbf{职称：}}
\newcommand\zhicheng[1]{\def\zhchstring{#1}}
\newcommand\jxrq{\textbf{教学日期：}}
\newcommand\riqi[1]{\def\rqstring{#1}}
\newcommand\dx{\textbf{授课对象：}}
\newcommand\duixiang[1]{\def\dxstring{#1}}
\newcommand\rsh{\textbf{实验人数：}}
\newcommand\renshu[1]{\def\rshstring{#1}}
\newcommand\lx{\textbf{实验类型（验证型、综合型、设计型、创新型）：}}
\newcommand\leixing[1]{\def\lxstring{#1}}
\newcommand\fz{\textbf{实验分组：}}
\newcommand\fenzu[1]{\def\fzstring{#1}}
\newcommand\xsh{\textbf{学时数：}}
\newcommand\xueshi[1]{\def\xshstring{#1}}
\newcommand\jc{\textbf{教材版本：}}
\newcommand\jiaocai[1]{\def\jcstring{#1}}


%每一页头部类似于页眉的部分
%首页
\newcommand\firstHeader{
	\begin{center}
		\textbf{\Large 天\ 津\ 医\ 科\ 大\ 学\ \uline{实\ 验\ 课}\ 教\ 案\ 首\ 页}
	\end{center}
	\begin{flushright}
	\small（共 \pageref{LastPage} 页、第 \thepage 页）
	\end{flushright}
	\vspace*{-10pt}
	\hrule
}
%续页
\newcommand\otherHeader{
	\begin{center}
		\textbf{\Large 天\ 津\ 医\ 科\ 大\ 学\ \uline{实\ 验\ 课}\ 教\ 案\ 续\ 页}
	\end{center}
	\begin{flushright}
	\small（共 \pageref{LastPage} 页、第 \thepage 页）
	\end{flushright}
	\vspace*{-10pt}
	\hrule
	\vspace*{5pt}
}

%首页底部类似于页脚的部分
%首页
\newcommand\firstTail{
	\vfill
	\smallskip
	\hrule
	\smallskip
	{\bf \noindent 主任签字：\hfill \quad 年 \quad 月 \quad 日 \hfill 教务处制}
}
%续页
\newcommand\otherTail{
	\vfill
	\hrule
	\vspace*{-10pt}
	\begin{flushright}
	\textbf{教务处制}
	\end{flushright}
}


%首页上的课程信息部分
\renewcommand\maketitle{
	\smallskip
	\noindent \kch \kchstring \hfill \shy \shystring
	\smallskip
	\hrule
	\smallskip
	\noindent \jsh \jshstring \hfill \zhch \zhchstring \hfill \jxrq \rqstring
	\smallskip
	\hrule
	\smallskip
	\noindent \dx \dxstring \hfill \rsh \rshstring
	\smallskip
	\hrule
	\smallskip
	\noindent \lx \lxstring \hfill \fz \fzstring
	\smallskip
	\hrule
	\smallskip
	\noindent \xsh \xshstring \hfill \jc \jcstring
}
%调整表格内行距
%\renewcommand\arraystretch{1.1}
\newcommand\maketitleTable{
\begin{tabu} to 0.752\textwidth {X[1.6,l]X[1,l]X[3,l]}
	\noalign{\smallskip}
	\multicolumn{2}{l}{\kch \kchstring} & \shy \shystring\\
	\hline
	\noalign{\smallskip}
	\jsh \jshstring & \zhch \zhchstring & \jxrq \rqstring\\
	\hline
	\noalign{\smallskip}
	\multicolumn{2}{l}{\dx \dxstring} & \rsh \rshstring\\
	\hline
	\noalign{\smallskip}
	\multicolumn{2}{l}{\lx \lxstring} & \fz \fzstring\\
	\hline
	\noalign{\smallskip}
	\xsh \xshstring & \multicolumn{2}{l}{\jc \jcstring}\\
\end{tabu}
}


%首页上的各个模块
\newenvironment{module}[1]{
	\smallskip
	\hrule
	\smallskip
	\noindent
}{\smallskip}

%实验目的与要求
\newcommand\mudi[1]{
\begin{module}
	\noindent \textbf{实验目的与要求：}
	{#1}
\end{module}
}
%实验内容及学时分配
\newcommand\fenpei[1]{
\begin{module}
	\noindent \textbf{实验内容及学时分配：}
	{#1}
\end{module}
}
%仪器和材料
\newcommand\cailiao[1]{
\begin{module}
	\noindent \textbf{主要仪器和实验材料：}
	{#1}
\end{module}
}
%重点和难点
\newcommand\zhongdian[1]{
\begin{module}
	\noindent \textbf{实验重点、难点及解决策略：}
	{#1}
\end{module}
}
%思考题
\newcommand\sikao[1]{
\begin{module}
	\noindent \textbf{思考题：}
	{#1}
\end{module}
}
%参考资料
\newcommand\cankao[1]{
\begin{module}
	\noindent \textbf{参考资料：}
	{#1}
\end{module}
}


