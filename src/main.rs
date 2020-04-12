#![windows_subsystem = "windows"]

extern crate meval;

use iced::{
    button, text_input, Align, Button, Column,
    Container, Element, Length, Row, Sandbox,
    Settings, Text, TextInput,
};
use std::io::ErrorKind;
use std::convert::TryFrom;

pub fn main() {
    Styling::run(Settings::default())
}

#[derive(Default)]
struct Styling {
    theme: style::Theme,

    function_input: text_input::State,
    function_input_value: String,

    steps_input: text_input::State,
    steps_input_value: String,

    a_input: text_input::State,
    a_input_value: String,

    b_input: text_input::State,
    b_input_value: String,

    x0_input: text_input::State,
    x0_input_value: String,

    y0_input: text_input::State,
    y0_input_value: String,

    result_text_value: String,

    euler_button: button::State,
    runge_button: button::State,
    adams_button: button::State,
}

#[derive(Debug, Clone)]
enum Message {
    FunctionInputChanged(String),

    StepsInputChanged(String),
    AInputChanged(String),
    BInputChanged(String),
    X0InputChanged(String),
    Y0InputChanged(String),

    EulerButtonPressed,
    RungeButtonPressed,
    AdamsButtonPressed,
}

struct Parsed {
    steps: i32,

    x0: f64,
    y0: f64,

    start: f64,
    end: f64,

    expr: meval::Expr,
}

impl TryFrom<&Styling> for Parsed {
    type Error = ErrorKind;

    fn try_from(source: &Styling) -> Result<Self, Self::Error> {
        let expr = source.function_input_value.parse();
        if expr.is_err() {
            return Err(ErrorKind::InvalidData);
        }

        Ok(Parsed {
            steps: source.steps_input_value.parse().unwrap_or(0i32),
            x0: source.x0_input_value.parse().unwrap_or(0f64),
            y0: source.y0_input_value.parse().unwrap_or(0f64),
            start: source.a_input_value.parse().unwrap_or(0f64),
            end: source.b_input_value.parse().unwrap_or(0f64),
            expr: expr.unwrap(),
        })
    }
}

fn format_result(title: String, iterations: &Vec<Vec<f64>>, steps: i32, h: f64) -> String {
    let mut iterations_string = String::from(format!("Method {}\n", title));
    iterations_string.push_str(format!("Steps: {} (h = {})\n", steps, h).as_str());
    iterations_string.push_str("#    x           y\n");

    let mut i = 0;
    for ar in iterations {
        let x = ar[0];
        let y = ar[1];

        iterations_string.push_str(format!("{}    {:.2}      {:.6}\n", i, x, y).as_str());

        i += 1;
    }

    iterations_string
}

fn fix_xd(inp: f64) -> f64 {
    (f64::round(inp * 10000f64) as i64) as f64 / 10000f64
}

fn euler(parsed: Parsed) -> String {
    // Converting data to numbers from strings
    let func = parsed.expr.bind2("x", "y").unwrap();

    // Algorithm start
    let mut iterations: Vec<Vec<f64>> = Vec::new();

    let mut xy = Vec::new();
    xy.push(parsed.x0);
    xy.push(parsed.y0);
    iterations.push(xy);

    let h: f64 = (parsed.end - parsed.start) / parsed.steps as f64;

    for i in 0..parsed.steps {
        let mut iter = Vec::new();

        let p_itr = iterations.get(i as usize).unwrap(); //previous iteration
        let px = *p_itr.get(0).unwrap();
        let py = *p_itr.get(1).unwrap();
        let pf = func(px, py);

        iter.push(fix_xd(px + h));
        iter.push(py + h / 2f64 * (pf + func(px, py + h * pf)));

        iterations.push(iter);
    }

    format_result(String::from("Euler"), &iterations, parsed.steps, h)
}

fn runge(parsed: Parsed) -> String {
    // Converting data to numbers from strings
    let func = parsed.expr.bind2("x", "y").unwrap();

    // Algorithm start
    let mut iterations: Vec<Vec<f64>> = Vec::new();

    let mut xy = Vec::new();
    xy.push(parsed.x0);
    xy.push(parsed.y0);
    iterations.push(xy);

    let h: f64 = (parsed.end - parsed.start) / parsed.steps as f64;

    for i in 0..parsed.steps {
        let mut iter = Vec::new();
        let p_iter = iterations.get(i as usize).unwrap(); //previous iteration

        let px = *p_iter.get(0).unwrap();
        let py = *p_iter.get(1).unwrap();

        let k1 = func(px, py);
        let k2 = func(px + h / 2f64, py + (h * k1) / 2f64);
        let k3 = func(px + h / 2f64, py + (h * k2) / 2f64);
        let k4 = func(px + h, py + (h * k3));

        iter.push(fix_xd(px + h));
        iter.push(py + h * (k1 + (2f64 * k2) + (2f64 * k3) + k4) / 6f64);
        iterations.push(iter);
    }

    format_result(String::from("Runge-Kutta"), &iterations, parsed.steps, h)
}

fn adams(parsed: Parsed) -> String {
    // Converting data to numbers from strings
    let func = parsed.expr.bind2("x", "y").unwrap();

    // Algorithm start
    let mut iterations: Vec<Vec<f64>> = Vec::new();
    let h: f64 = (parsed.end - parsed.start) / parsed.steps as f64;

    let mut xy = Vec::new();
    xy.push(parsed.x0);
    xy.push(parsed.y0);
    iterations.push(xy);

    for i in 1..3 {
        let mut iter = Vec::new();
        let p_iter = iterations.get(i as usize - 1).unwrap(); //previous iteration

        let px = *p_iter.get(0).unwrap();
        let py = *p_iter.get(1).unwrap();

        iter.push(h * i as f64);
        iter.push(py + h * func(px, py));

        iterations.push(iter);
    }

    for i in 2..parsed.steps {
        let mut iter = Vec::new();

        let p1y = *iterations.get(i as usize).unwrap().get(1).unwrap();
        let p2y = *iterations.get(i as usize - 1).unwrap().get(1).unwrap();
        let p3y = *iterations.get(i as usize - 2).unwrap().get(1).unwrap();

        let t = parsed.start + i as f64 * h;
        let k1 = 23f64 * func(t, p1y);
        let k2 = 16f64 * func(t - h, p2y);
        let k3 = 5f64 * func(t - 2f64 * h, p3y);

        iter.push(t + h);
        iter.push(p1y + (h / 12f64) * (k1 - k2 + k3));

        iterations.push(iter);
    }

    format_result(String::from("Adams"), &iterations, parsed.steps, h)
}

impl Sandbox for Styling {
    type Message = Message;

    fn new() -> Self {
        Styling::default()
    }

    fn title(&self) -> String {
        String::from("Numerical Methods: Laboratory work №4")
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::FunctionInputChanged(value) => self.function_input_value = value,
            Message::StepsInputChanged(value) => self.steps_input_value = value,
            Message::AInputChanged(value) => self.a_input_value = value,
            Message::BInputChanged(value) => self.b_input_value = value,
            Message::X0InputChanged(value) => self.x0_input_value = value,
            Message::Y0InputChanged(value) => self.y0_input_value = value,

            Message::EulerButtonPressed => {
                let parsed = Parsed::try_from(& *self);
                if parsed.is_err() {
                    &self.result_text_value.clear();
                    &self.result_text_value.push_str("An error occurred.");
                    return;
                }

                let result = euler(parsed.unwrap());
                &self.result_text_value.clear();
                &self.result_text_value.push_str(result.as_str());
            }
            Message::RungeButtonPressed => {
                let parsed = Parsed::try_from(& *self);
                if parsed.is_err() {
                    &self.result_text_value.clear();
                    &self.result_text_value.push_str("An error occurred.");
                    return;
                }

                let result = runge(parsed.unwrap());
                &self.result_text_value.clear();
                &self.result_text_value.push_str(result.as_str());
            }
            Message::AdamsButtonPressed => {
                let parsed = Parsed::try_from(& *self);
                if parsed.is_err() {
                    &self.result_text_value.clear();
                    &self.result_text_value.push_str("An error occurred.");
                    return;
                }

                let result = adams(parsed.unwrap());
                &self.result_text_value.clear();
                &self.result_text_value.push_str(result.as_str());
            }
        }
    }

    fn view(&mut self) -> Element<Message> {
        let function_input = TextInput::new(
            &mut self.function_input,
            "Enter your function...",
            &self.function_input_value,
            Message::FunctionInputChanged,
        ).padding(10).size(20).style(self.theme);

        let steps_input = TextInput::new(
            &mut self.steps_input,
            "Number of steps",
            &self.steps_input_value,
            Message::StepsInputChanged,
        ).padding(10).size(20).style(self.theme);

        let a_input = TextInput::new(
            &mut self.a_input,
            "Start",
            &self.a_input_value,
            Message::AInputChanged,
        ).padding(10).size(20).style(self.theme);

        let b_input = TextInput::new(
            &mut self.b_input,
            "End",
            &self.b_input_value,
            Message::BInputChanged,
        ).padding(10).size(20).style(self.theme);

        let x0_input = TextInput::new(
            &mut self.x0_input,
            "x",
            &self.x0_input_value,
            Message::X0InputChanged,
        ).padding(10).size(20).style(self.theme);

        let y0_input = TextInput::new(
            &mut self.y0_input,
            "y",
            &self.y0_input_value,
            Message::Y0InputChanged,
        ).padding(10).size(20).style(self.theme);

        let euler_button = Button::new(
            &mut self.euler_button,
            Text::new("Euler Method"),
        ).padding(10).on_press(Message::EulerButtonPressed).style(self.theme);

        let runge_button = Button::new(
            &mut self.runge_button,
            Text::new("Runge–Kutta Method"),
        ).padding(10).on_press(Message::RungeButtonPressed).style(self.theme);

        let adams_button = Button::new(
            &mut self.adams_button,
            Text::new("Adams Method"),
        ).padding(10).on_press(Message::AdamsButtonPressed).style(self.theme);

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .max_width(600)
            .align_items(Align::Center)
            .push(
                Column::new().push(Row::new().spacing(10).align_items(Align::Center)
                    .push(Text::new("Function"))
                    .push(function_input)
                )
                    .push(Row::new().spacing(10).align_items(Align::Center)
                        .push(Text::new("Steps"))
                        .push(steps_input)
                    )
                    .push(Row::new().spacing(10).align_items(Align::Center)
                        .push(Text::new("a"))
                        .push(a_input)
                        .push(Text::new("b"))
                        .push(b_input)
                    )
                    .push(Row::new().spacing(10).align_items(Align::Center)
                        .push(Text::new("x0"))
                        .push(x0_input)
                        .push(Text::new("y0"))
                        .push(y0_input)
                    )
                    .push(Row::new().spacing(10).align_items(Align::Center)
                        .push(euler_button)
                        .push(runge_button)
                        .push(adams_button)
                    )
                    .push(Row::new().spacing(10).align_items(Align::Center)))
            .push(Column::new().push(Text::new(&self.result_text_value)));

        Container::new(content)
            .width(Length::Fill).height(Length::Fill)
            .center_x().center_y()
            .style(self.theme)
            .into()
    }
}

mod style {
    use iced::{button, checkbox, container, progress_bar, radio, scrollable, slider, text_input};

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Theme {
        Light
    }

    impl Theme {
        pub const ALL: [Theme; 1] = [Theme::Light];
    }

    impl Default for Theme {
        fn default() -> Theme {
            Theme::Light
        }
    }

    impl From<Theme> for Box<dyn container::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    impl From<Theme> for Box<dyn radio::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    impl From<Theme> for Box<dyn text_input::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    impl From<Theme> for Box<dyn button::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => light::Button.into()
            }
        }
    }

    impl From<Theme> for Box<dyn scrollable::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    impl From<Theme> for Box<dyn slider::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    impl From<Theme> for Box<dyn progress_bar::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    impl From<Theme> for Box<dyn checkbox::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default()
            }
        }
    }

    mod light {
        use iced::{button, Background, Color, Vector};

        pub struct Button;

        impl button::StyleSheet for Button {
            fn active(&self) -> button::Style {
                button::Style {
                    background: Some(Background::Color(Color::from_rgb(
                        0.11, 0.42, 0.87,
                    ))),
                    border_radius: 12,
                    shadow_offset: Vector::new(1.0, 1.0),
                    text_color: Color::from_rgb8(0xEE, 0xEE, 0xEE),
                    ..button::Style::default()
                }
            }

            fn hovered(&self) -> button::Style {
                button::Style {
                    text_color: Color::WHITE,
                    shadow_offset: Vector::new(1.0, 2.0),
                    ..self.active()
                }
            }
        }
    }
}