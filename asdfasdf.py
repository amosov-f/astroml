from docxtpl import DocxTemplate


def main():
    doc = DocxTemplate("Заявление_о_присоединении.docx")
    context = {'myvariable': "World company"}
    doc.render(context)
    doc.save("generated_doc.docx")

if __name__ == '__main__':
    main()