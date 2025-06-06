#import "@local/gost732-2017:0.3.0": *
#import "@local/bmstu:0.2.0": *

#show: приложение.with(буква: "А", содержание: [ Техническое задание ])

// TODO Поменять количество страниц в прил А на титульнике
#страница(image("../титулы/тз-титул.jpg", height: 100%), номер: нет)

#ненумерованный_заголовок(содержание: нет)[ 1 ВВЕДЕНИЕ ]
#v(1em)

Настоящее техническое задание распространяется на разработку цифрового омметра, именуемого в дальнейшем «устройством». Данное устройство принимает подключаемый измеряемый резистор и отображает его электрическое сопротивление на цифровом индикаторе. Устройство должно быть выполнено на элементной базе КМОП.


#ненумерованный_заголовок(содержание: нет)[ 2 ОСНОВАНИЯ ДЛЯ РАЗРАБОТКИ ]
#v(1em)

Устройство разрабатывается в качестве курсового проекта на основе
плана учебной работы студентов МГТУ им. Баумана кафедры ИУ6 «Компьютерные системы и сети».

#ненумерованный_заголовок(содержание: нет)[ 3 НАЗНАЧЕНИЕ РАЗРАБОТКИ ]
#v(1em)

Устройство предназначено для измерения электрического сопротивления подключаемых резисторов в заданном диапазоне и отображения результата измерения в цифровом виде. Измерение основано на определении времени заряда конденсатора известной емкости через измеряемое сопротивление до заданного порогового напряжения.


#ненумерованный_заголовок(содержание: нет)[ 4 ЦЕЛИ И РЕШАЕМЫЕ ЗАДАЧИ ]
#v(1em)

4.1 Цель работы

Целью курсового проектирования является разработка цифрового омметра на основе измерения времени заряда RC-цепи с использованием элементной базы КМОП.

4.2 Решаемые задачи

4.2.1 Анализ технического задания и возможных путей решения
поставленной задачи.

4.2.2 Обоснование и синтез электрической функциональной схемы устройства.

4.2.3 Выбор элементной базы на основании технических требований. 

4.2.4 Разработка электрической принципиальной схемы устройства. 

4.2.5 Построение временных диаграмм.

4.2.6 Расчет параметров быстродействия и мощности устройства.

#ненумерованный_заголовок(содержание: нет)[ 5 ТЕХНИЧЕСКИЕ ТРЕБОВАНИЯ ]
#v(1em)

К разрабатываемому устройству предъявляются следующие
требования к составу и параметру технических средств:

5.1 Устройство должно измерять электрическое сопротивление в диапазоне от 1Ом до 1МОм. Измерение должно производиться методом определения времени заряда RC-цепи до порогового напряжения.

5.2 Логика элементов – КМОП.

5.3 Тактовая частота генератора импульсов – 1 МГц.

5.4  Мощность потребления – не более 10 Вт.

5.5 К разрабатываемому устройству предъявляются условия
эксплуатации в соответствие с СанПиН2.2.2/2.4.1340-03.
Требования к маркировке и упаковке разрабатываемого устройства
не предъявляются.

5.6 Требования к транспортированию и хранению разрабатываемого устройства не предъявляются. \ \

#ненумерованный_заголовок(содержание: нет)[ 6 ТРЕБОВАНИЯ К ДОКУМЕНТАЦИИ ]
#v(1em)

6.1 В состав сопровождающей документации должны входить:

6.1.1 Расчетно-пояснительная записка на 20 -- 30 листах формата А4

6.1.2 Техническое задание (Приложение А)

6.1.3 Спецификация (Приложение Е)

6.2 Графическая часть должна быть включена в расчетно-пояснительную записку в качестве приложений и иллюстраций:

6.2.1 Временные диаграммы (Приложение Г).

6.2.2 Схема структурная (Приложение Б).

6.2.3 Схема электрическая функциональная (Приложение B).

6.2.4 Схема электрическая принципиальная (Приложение Д).

#ненумерованный_заголовок(содержание: нет)[ 7 СТАДИИ И ЭТАПЫ РАЗРАБОТКИ ]
#v(1em)

Этапы разработки курсовой работы указаны в таблице @этапы-разработки.

#таблица(table(
  columns: (auto, auto, auto, auto),
  [ № ], [ Название этапа ], [ Срок, % выполнения ], [ Отчетность ],
  [ 1 ], [ Выдача задания на проект ], [ - ], [ - ],
  [ 2 ], [ Анализ технического задания. Обоснование и синтез электрической функциональной схемы узла. ], [ 3 нед., 20% ], [ Функциональная схема ],
  [ 3 ], [ Обоснование выбора элементной базы и разработка электрической принципиальной схемы узла. ], [ 6 нед., 40% ], [ Принципиальная схема ],
  [ 4 ], [ Разработка временных диаграмм функционирования узла. Выполнение расчетов. ], [ - ], [ Временные диаграммы. Расчеты ],
  [ 5 ], [ Смотр состояния проекта. ], [ - ], [ - ],
  [ 6 ], [ Конструкторское проектирование печатной платы. Выполнение конструкторских расчетов. ], [ 11 нед., 90% ], [ - ],
  [ 7 ], [ Окончательное оформление графической части проекта и расчетно-пояснительной записки. ], [ 12 нед., 100% ], [ Расчетно-поясни- тельная записка. ],
  [ 8 ], [ Подготовка доклада и защита курсового проекта. ], [ 13 нед., - ], [ Доклад (3-5 минут). ],
))[ Этапы разработки ] <этапы-разработки>

#ненумерованный_заголовок(содержание: нет)[ 8 ПОРЯДОК КОНТРОЛЯ И ПРИЕМА ]
#v(1em)

8.1 Порядок контроля

Контроль выполнения осуществляется руководителем еженедельно. 

8.2 Порядок защиты

Защита осуществляется перед комиссией, состоящей из преподавателей кафедры ИУ6.

8.3 Срок защиты

Срок защиты: 15-16 недели.

#ненумерованный_заголовок(содержание: нет)[ 9 ПРИМЕЧАНИЕ ]
#v(1em)

В процессе выполнения работы возможно уточнение отдельных
требований технического задания по взаимному согласованию руководителя и исполнителя.