#import "@local/gost732-2017:0.3.0": *
#import "@local/bmstu:0.1.1": *

#show: гост732-2017


#страница(image("титулы/титул-рпз.jpg", height: 100%), номер: нет)

#страница(image("титулы/бланк-задание.jpg", height: 100%), номер: нет)

#include "разделы/0-реферат.typ"


#include "разделы/1-сокращения.typ"

#содержание()

#include "разделы/2-введение.typ"

#include "разделы/3-анализ-требований.typ"

#include "разделы/4-проектрирование-функ.typ"

#include "разделы/5-построение-временных.typ"

#include "разделы/6-проектрирование-принцип.typ"

#include "разделы/7-быстродействие.typ"

#include "разделы/8-мощность.typ"

#include "разделы/9-заключение.typ"


#bibliography("liba.yml")


//#include "разделы/прил-а.typ"
// 
//#include "разделы/прил-Б.typ"

//#include "разделы/прил-В.typ"

//#include "разделы/прил-Г.typ"

//#include "разделы/прил-Д.typ"

//#include "разделы/прил-Е.typ"
