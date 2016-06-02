# 该 Makefile 用于编译一个 Beamer 制作的幻灯片
# 运行环境：Linux、texlive2010、xpdf 阅读器

PDFC = xelatex
READER = evince

TARGET = car.pdf
OUTGROWTH = *.out *.log *.aux *.nav *.snm *.toc
SOURCE = car.tex

$(TARGET): $(SOURCE)
	$(PDFC)	$(SOURCE)
	$(PDFC)	$(SOURCE)

read:
	$(READER) $(TARGET)

clean:
	-rm $(OUTGROWTH)

cleanall:
	-rm $(OUTGROWTH) $(TARGET)

.PHONY:read clean cleanall
